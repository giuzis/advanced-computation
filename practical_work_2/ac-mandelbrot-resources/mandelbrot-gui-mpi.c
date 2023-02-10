#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <mpi.h>

////////////////////////////////////////////////////////////////////////
//user defined datatypes
typedef struct {unsigned char r, g, b;} rgb_t;

typedef struct {
	int quit, color_rotate, height, invert, max_iter, refresh, saturation, tex_h, tex_size, tex_w, width;
	double cx, cy, scale;
} global_parameters_t;

MPI_Datatype PARAMETERS_TYPE; // global_parameters_t

////////////////////////////////////////////////////////////////////////
//global variables
//(do not change from the defaults except if solicited in the work text)
global_parameters_t GLOBAL_parameters = {
	   0,  // int GLOBAL_quit = 0;
	   0,  // int GLOBAL_color_rotate = 0;
	   0,  // int GLOBAL_height;
	   0,  // int GLOBAL_invert = 0;
	 256,  // int GLOBAL_max_iter = 256;
	   1,  // int GLOBAL_refresh = 1;
	   1,  // int GLOBAL_saturation = 1;
	   0,  // int GLOBAL_tex_h;
	   0,  // int GLOBAL_tex_size = 0;
	   0,  // int GLOBAL_tex_w;
	   0,  // int GLOBAL_width;
	-0.6,  // double GLOBAL_cx = -.6;
	   0,  // double GLOBAL_cy = 0;
  1./256   // double GLOBAL_scale = 1;
};

rgb_t **GLOBAL_tex = 0;

int GLOBAL_gwin;
GLuint GLOBAL_texture;

// GLOBAL_zoomin contains the path (sequence of <x,y> mouse coordinates) to zoom in the fractal)
#define GLOBAL_zoomin_num_pairs 36
int GLOBAL_zoomin[GLOBAL_zoomin_num_pairs]={538,237,491,369,522,383,492,372,504,369,510,385,491,353,472,329,516,392,518,393,478,327,506,400,501,320,487,361,501,363,501,363,486,363,486,363};

int GLOBAL_window_width=1024;
int GLOBAL_window_height=768;

// GLOBAL mpi variables
int GLOBAL_numtasks, GLOBAL_rank;

////////////////////////////////////////////////////////////////////////
//function prototypes
void render();
void keypress(unsigned char key, int x, int y);
void mouseclick(int button, int state, int x, int y);
void resize(int w, int h);
void init_gfx(int *c, char **v);

void alloc_tex();
int set_texture();

void hsv_to_rgb(int hue, int min, int max, rgb_t *p);
void calc_mandel();

void print_menu();
void screen_dump();

////////////////////////////////////////////////////////////////////////
// function implementations

void print_menu()
{
	printf("\n\nkeys:\n\t"
	"q: quit	\n\t"
	"ESC: reset to initial frame\n\t"
	"r: color rotation\n\t"
	"c: monochrome\n\t"
	"s: screen dump\n\t"
	"<, >: decrease/increase max iteration\n\t"
	"I: max iteration=4096\n\t"
	"i: max iteration=128\n\t"
	"mouse buttons to zoom\n\t"
	"z: automatic zoom in (one step)\n\t"
	"Z: automatic zoom in (all steps)\n");
}

//////////////////////////////////////////////////////////////////////// 
void screen_dump()
{
	static int dump=1;
	char fn[100];
	int i;
	sprintf(fn, "screen%03d.ppm", dump++);
	FILE *fp = fopen(fn, "w");
	fprintf(fp, "P6\n%d %d\n255\n", GLOBAL_parameters.width, GLOBAL_parameters.height);
	for (i = GLOBAL_parameters.height - 1; i >= 0; i--)
		fwrite(GLOBAL_tex[i], 1, GLOBAL_parameters.width * 3, fp);
	fclose(fp);
	printf("%s written\n", fn);
}

////////////////////////////////////////////////////////////////////////
void hsv_to_rgb(int hue, int min, int max, rgb_t *p)
{
	if (min == max) max = min + 1;
	if (GLOBAL_parameters.invert) hue = max - (hue - min);
	if (!GLOBAL_parameters.saturation) {
		p->r = p->g = p->b = 255 * (max - hue) / (max - min);
		return;
	}
	double h = fmod(GLOBAL_parameters.color_rotate + 1e-4 + 4.0 * (hue - min) / (max - min), 6);
#	define VAL 255
	double c = VAL * GLOBAL_parameters.saturation;
	double X = c * (1 - fabs(fmod(h, 2) - 1));
 
	p->r = p->g = p->b = 0;
 
	switch((int)h) {
	case 0: p->r = c; p->g = X; return;
	case 1:	p->r = X; p->g = c; return;
	case 2: p->g = c; p->b = X; return;
	case 3: p->g = X; p->b = c; return;
	case 4: p->r = X; p->b = c; return;
	default:p->r = c; p->b = X;
	}
}

////////////////////////////////////////////////////////////////////////
void calc_mandel(){	
	int min = GLOBAL_parameters.max_iter, max = 0, i, j;
	long double x, y, zx, zy, zx2, zy2;
	rgb_t *px;

	int height = (GLOBAL_parameters.height / GLOBAL_numtasks);
	int rest = (GLOBAL_parameters.height % GLOBAL_numtasks);

	for (i = 0; i < (height + (rest * (!GLOBAL_rank))); i++) {
		px = GLOBAL_tex[i];
		y = ((GLOBAL_rank * height) + (rest * (GLOBAL_rank != 0)) + i - GLOBAL_parameters.height / 2) * GLOBAL_parameters.scale + GLOBAL_parameters.cy;
		for (j = 0; j  < GLOBAL_parameters.width; j++, px++) {
			x = (j - GLOBAL_parameters.width / 2) * GLOBAL_parameters.scale + GLOBAL_parameters.cx;
			int iter = 0;
 
			zx = hypot(x - .25, y);
			if (x < zx - 2 * zx * zx + .25) iter = GLOBAL_parameters.max_iter;
			if ((x + 1)*(x + 1) + y * y < 1/16) iter = GLOBAL_parameters.max_iter;
 
			zx = zy = zx2 = zy2 = 0;
			for (; iter < GLOBAL_parameters.max_iter && zx2 + zy2 < 4; iter++) {
				zy = 2 * zx * zy + y;
				zx = zx2 - zy2 + x;
				zx2 = zx * zx;
				zy2 = zy * zy;
			}
			if (iter < min) min = iter;
			if (iter > max) max = iter;
			*(unsigned short*)px = iter;
		}
	}

	MPI_Allreduce(&min, &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&max, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	for (i = 0; i < (height + (rest * (!GLOBAL_rank))); i++){
		px = GLOBAL_tex[i];
		for (j = 0; j  < GLOBAL_parameters.width; j++, px++)
			hsv_to_rgb(*(unsigned short*)px, min, max, px);			
	}

	int size = height * GLOBAL_parameters.width * 3;

	if (GLOBAL_rank == 0)
		MPI_Gather(MPI_IN_PLACE, -1, NULL, GLOBAL_tex[rest], size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	else
		MPI_Gather(GLOBAL_tex[0], size, MPI_UNSIGNED_CHAR, NULL, -1, NULL, 0, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////
void alloc_tex(){
	int i, ow = GLOBAL_parameters.tex_w, oh = GLOBAL_parameters.tex_h; //backup current texture dimensions
	int tex_w = 1, tex_h = 1;

	int new_height = GLOBAL_parameters.height;

	if (GLOBAL_rank != 0) new_height /= GLOBAL_numtasks; // get the height for each task that is not rank 0

	while(tex_w < GLOBAL_parameters.width) tex_w <<= 1; // smallest power of two >= width
	while(tex_h < new_height) tex_h <<= 1; // smallest power of two >= height

	GLOBAL_parameters.tex_w = tex_w;
	GLOBAL_parameters.tex_h = tex_h;

	if (GLOBAL_parameters.tex_h != oh || GLOBAL_parameters.tex_w != ow) { // if the dimensions are the different than before, realocate
		GLOBAL_parameters.tex_size = GLOBAL_parameters.tex_h * sizeof(rgb_t*) + GLOBAL_parameters.tex_h * GLOBAL_parameters.tex_w * 3;
		GLOBAL_tex = realloc(GLOBAL_tex, GLOBAL_parameters.tex_size);

		for (GLOBAL_tex[0] = (rgb_t *)(GLOBAL_tex + GLOBAL_parameters.tex_h), i = 1; i < GLOBAL_parameters.tex_h; i++)
			GLOBAL_tex[i] = GLOBAL_tex[i - 1] + GLOBAL_parameters.tex_w;   // uses rgb_t* arithmetic pointer, where each unit corresponds to 3 bytes
	}
}

////////////////////////////////////////////////////////////////////////~
int set_texture(){
	MPI_Bcast(&GLOBAL_parameters, 1, PARAMETERS_TYPE, 0, MPI_COMM_WORLD);
	
	if (GLOBAL_parameters.quit){
		
		return 0;
	}
	if (GLOBAL_parameters.refresh) alloc_tex();
	
	calc_mandel();
 
	if (GLOBAL_parameters.refresh && GLOBAL_rank == 0) {
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, GLOBAL_texture);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, GLOBAL_parameters.tex_w, GLOBAL_parameters.tex_h, 0, GL_RGB, GL_UNSIGNED_BYTE, GLOBAL_tex[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		render();
	}

	return 1;
}

////////////////////////////////////////////////////////////////////////
void resize(int w, int h)
{
	GLOBAL_parameters.width = w;
	GLOBAL_parameters.height = h;
 
	glViewport(0, 0, w, h);
	glOrtho(0, w, 0, h, -1, 1);
 
	set_texture();
}

////////////////////////////////////////////////////////////////////////
void mouseclick(int button, int state, int x, int y)
{
	if (state != GLUT_UP) return;
	
	GLOBAL_parameters.cx += (x - GLOBAL_parameters.width / 2) * GLOBAL_parameters.scale;
	GLOBAL_parameters.cy -= (y - GLOBAL_parameters.height/ 2) * GLOBAL_parameters.scale;
 
	switch(button) {
	case GLUT_LEFT_BUTTON: /* zoom in */
		if (GLOBAL_parameters.scale > fabs(x) * 1e-16 && GLOBAL_parameters.scale > fabs(y) * 1e-16)
			GLOBAL_parameters.scale /= 2;
		break;
	case GLUT_RIGHT_BUTTON: /* zoom out */
		GLOBAL_parameters.scale *= 2;
		break;
	/* any other button recenters */
	}
	set_texture();
	//print_menu(); // uncomment for convenience; comment for benchmarking
}

////////////////////////////////////////////////////////////////////////
void keypress(unsigned char key, int x, int y){
	static int zoomin_x=0, zoomin_y=1; // RUF
	// Ruf: where to start fetching mouse coordinates from GLOBAL_zoomin
	//      (first x coordinate is at GLOBAL_zoomin[0])
	//      (first y coordinate is at GLOBAL_zoomin[1])
	//      (next coordinates are at distance 2: see +=2 at 'z' and 'Z" bellow)
	
	switch(key) {
	case 'q': 
            //  glutDestroyWindow(GLOBAL_gwin);
            //  glFinish();
            //  free(GLOBAL_tex);
			 GLOBAL_parameters.quit = 1;
			 MPI_Bcast(&GLOBAL_parameters, 1, PARAMETERS_TYPE, 0, MPI_COMM_WORLD);
			 glutLeaveMainLoop();
			 break;
	
	case 27: // Esc
	         GLOBAL_parameters.scale = 1./256;
	         GLOBAL_parameters.cx = -.6;
	         GLOBAL_parameters.cy = 0; 
	         break;
 
	case 'r':
             GLOBAL_parameters.color_rotate = (GLOBAL_parameters.color_rotate + 1) % 6;
              break;

	case '>':
    case '.':
             GLOBAL_parameters.max_iter += 128;
             if (GLOBAL_parameters.max_iter > 1 << 15) GLOBAL_parameters.max_iter = 1 << 15;
             printf("max iter: %d\n", GLOBAL_parameters.max_iter);
             break;
 
	case '<':
    case ',':
             GLOBAL_parameters.max_iter -= 128;
             if (GLOBAL_parameters.max_iter < 128) GLOBAL_parameters.max_iter = 128;
             printf("max iter: %d\n", GLOBAL_parameters.max_iter);
             break;

	case 'c':
             GLOBAL_parameters.saturation = 1 - GLOBAL_parameters.saturation;
             break;
 
	case 's':screen_dump(); return;

	case 'I':GLOBAL_parameters.max_iter = 4096; break;
	
	case 'i':GLOBAL_parameters.max_iter = 256; break;
	
	case ' ':GLOBAL_parameters.invert = !GLOBAL_parameters.invert; break;
	
	case 'z':// simulate one mouse click in order to dive one time in zoomin
             GLOBAL_parameters.refresh=1;
             mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]);
             zoomin_x+=2; zoomin_y+=2;
             break;

	case 'Z':// simulate many mouse clicks in order to dive fully in zoomin
             GLOBAL_parameters.refresh=1; // use 0 to avoid refreshing all but the last one
             for (zoomin_x=0, zoomin_y=1; zoomin_x < GLOBAL_zoomin_num_pairs; zoomin_x+=2, zoomin_y +=2) {
                 if (zoomin_x == GLOBAL_zoomin_num_pairs-2) GLOBAL_parameters.refresh=1;
                 mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]); 					
             }

             // simulate case 's'
             keypress('s', -1, -1);
              
             // simulate case 'q'
             keypress('q', -1, -1);
             return;
	}
	if (!GLOBAL_parameters.quit){
		set_texture();
		print_menu();
	}
}

////////////////////////////////////////////////////////////////////////
void render()
{
	double	x = (double)GLOBAL_parameters.width /GLOBAL_parameters.tex_w,
		    y = (double)GLOBAL_parameters.height/GLOBAL_parameters.tex_h;
 
	glClear(GL_COLOR_BUFFER_BIT);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
	glBindTexture(GL_TEXTURE_2D, GLOBAL_texture);
 
	glBegin(GL_QUADS);
 
	glTexCoord2f(0, 0); glVertex2i(0, 0);
	glTexCoord2f(x, 0); glVertex2i(GLOBAL_parameters.width, 0);
	glTexCoord2f(x, y); glVertex2i(GLOBAL_parameters.width, GLOBAL_parameters.height);
	glTexCoord2f(0, y); glVertex2i(0, GLOBAL_parameters.height);
 
	glEnd();
 
	glFlush();
	glFinish();
}

////////////////////////////////////////////////////////////////////////

void init_gfx(int *c, char **v){
	glutInit(c, v);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(GLOBAL_window_width, GLOBAL_window_height);
    
	GLOBAL_gwin = glutCreateWindow("Mandelbrot");
	glutDisplayFunc(render);
 
	glutKeyboardFunc(keypress);
	glutMouseFunc(mouseclick);
	glutReshapeFunc(resize);
	glGenTextures(1, &GLOBAL_texture);
	glutSetOption( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS); 
	set_texture();
}

////////////////////////////////////////////////////////////////////////
void setup_parameter_type(){
    int blockcounts[2]; MPI_Aint offsets[2]; MPI_Datatype oldtypes[2];
	// setup blockcounts and oldtypes 
	blockcounts[0] = 11; oldtypes[0] = MPI_INT; 
	blockcounts[1] = 3; oldtypes[1] = MPI_DOUBLE; 
	
	// setup displacements
	global_parameters_t parameters; MPI_Aint base_address;
	MPI_Get_address(&parameters, &base_address); 
	MPI_Get_address(&parameters.quit, &offsets[0]); 
	MPI_Get_address(&parameters.cx, &offsets[1]);

	offsets[0] = MPI_Aint_diff(offsets[0], base_address);
	offsets[1] = MPI_Aint_diff(offsets[1], base_address);

	// create structured derived data type
	MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &PARAMETERS_TYPE);
	MPI_Type_commit(&PARAMETERS_TYPE);
}
////////////////////////////////////////////////////////////////////////
int main(int c, char **v){
	// init mpi
	MPI_Init(&c, &v);
	MPI_Comm_size(MPI_COMM_WORLD, &GLOBAL_numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &GLOBAL_rank);
	
	// init parameters
	setup_parameter_type();
	
	// create mpi datatype for struct
	if (GLOBAL_rank == 0){
		init_gfx(&c, v);
		print_menu();
		glutMainLoop();
	}
	else{
		while (set_texture());
	}
	printf("rank %d: done\n", GLOBAL_rank);
	free(GLOBAL_tex);
	MPI_Type_free(&PARAMETERS_TYPE);
	MPI_Finalize();
	return 0;
}
