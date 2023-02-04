#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

////////////////////////////////////////////////////////////////////////
//user defined datatypes
typedef struct {unsigned char r, g, b;} rgb_t;


////////////////////////////////////////////////////////////////////////
//global variables
//(do not change from the defaults except if solicited in the work text)

rgb_t **GLOBAL_tex = 0;
int GLOBAL_gwin;
GLuint GLOBAL_texture;
int GLOBAL_width, GLOBAL_height;
int GLOBAL_tex_w, GLOBAL_tex_h;
double GLOBAL_scale = 1./256;
double GLOBAL_cx = -.6, GLOBAL_cy = 0;
int GLOBAL_color_rotate = 0;
int GLOBAL_saturation = 1;
int GLOBAL_invert = 0;

// GLOBAL_zoomin contains the path (sequence of <x,y> mouse coordinates) to zoom in the fractal)
#define GLOBAL_zoomin_num_pairs 36
int GLOBAL_zoomin[GLOBAL_zoomin_num_pairs]={538,237,491,369,522,383,492,372,504,369,510,385,491,353,472,329,516,392,518,393,478,327,506,400,501,320,487,361,501,363,501,363,486,363,486,363};

int GLOBAL_window_width=1024;
int GLOBAL_window_height=768;
int GLOBAL_refresh=1;
int GLOBAL_max_iter = 256;
int GLOBAL_tex_size=0;

////////////////////////////////////////////////////////////////////////
//function prototypes
void render();
void keypress(unsigned char key, int x, int y);
void mouseclick(int button, int state, int x, int y);
void resize(int w, int h);
void init_gfx(int *c, char **v);

void alloc_tex();
void set_texture();

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
	fprintf(fp, "P6\n%d %d\n255\n", GLOBAL_width, GLOBAL_height);
	for (i = GLOBAL_height - 1; i >= 0; i--)
		fwrite(GLOBAL_tex[i], 1, GLOBAL_width * 3, fp);
	fclose(fp);
	printf("%s written\n", fn);
}

////////////////////////////////////////////////////////////////////////
void hsv_to_rgb(int hue, int min, int max, rgb_t *p)
{
	if (min == max) max = min + 1;
	if (GLOBAL_invert) hue = max - (hue - min);
	if (!GLOBAL_saturation) {
		p->r = p->g = p->b = 255 * (max - hue) / (max - min);
		return;
	}
	double h = fmod(GLOBAL_color_rotate + 1e-4 + 4.0 * (hue - min) / (max - min), 6);
#	define VAL 255
	double c = VAL * GLOBAL_saturation;
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
void calc_mandel() 
{	
	int i, j, iter, min, max;
	rgb_t *px;
	double x, y, zx, zy, zx2, zy2;
	min = GLOBAL_max_iter; max = 0;
	for (i = 0; i < GLOBAL_height; i++) {
		px = GLOBAL_tex[i];
		y = (i - GLOBAL_height/2) * GLOBAL_scale + GLOBAL_cy;
		for (j = 0; j  < GLOBAL_width; j++, px++) {
			x = (j - GLOBAL_width/2) * GLOBAL_scale + GLOBAL_cx;
			iter = 0;
 
			zx = hypot(x - .25, y);
			if (x < zx - 2 * zx * zx + .25) iter = GLOBAL_max_iter;
			if ((x + 1)*(x + 1) + y * y < 1/16) iter = GLOBAL_max_iter;
 
			zx = zy = zx2 = zy2 = 0;
			for (; iter < GLOBAL_max_iter && zx2 + zy2 < 4; iter++) {
				zy = 2 * zx * zy + y;
				zx = zx2 - zy2 + x;
				zx2 = zx * zx;
				zy2 = zy * zy;
			}
			if (iter < min) min = iter;
			if (iter > max) max = iter;
			*(unsigned short *)px = iter;
		}
	}
 
	for (i = 0; i < GLOBAL_height; i++)
		for (j = 0, px = GLOBAL_tex[i]; j  < GLOBAL_width; j++, px++)
			hsv_to_rgb(*(unsigned short*)px, min, max, px);			
}

////////////////////////////////////////////////////////////////////////
// this function is called when the window is resized, or when the texture is created
// it allocates the texture memory block (GLOBAL_tex) and sets the texture parameters
// the texture memory block is a 2D array of rgb_t values, but it is stored as a 1D array
// of pointers to the rows of the 2D array, followed by the 2D array itself
// this is done to make the texture memory block compatible with OpenGL's glTexImage2D()
// function, which requires a 1D array of pointers to the rows of the 2D array
// the texture memory block is allocated with a size that is a power of 2, so that
// the texture can be used with OpenGL's GL_REPEAT texture wrap mode
// the texture memory block is allocated with a size that is a multiple of 4, so that
// the texture can be used with OpenGL's GL_UNPACK_ALIGNMENT parameter
// the texture memory block is allocated with a size that is a multiple of 16, so that
// the texture can be used with OpenGL's GL_UNPACK_ROW_LENGTH parameter
// the texture memory block is allocated with a size that is a multiple of 64, so that
// the texture can be used with OpenGL's GL_UNPACK_SKIP_PIXELS and GL_UNPACK_SKIP_ROWS

void alloc_tex(){
	int i, ow = GLOBAL_tex_w, oh = GLOBAL_tex_h;
	       //backup current texture dimensions

    // if necessary, adjust texture dimensions to image dimensions 
	for (GLOBAL_tex_w = 1; GLOBAL_tex_w < GLOBAL_width;  GLOBAL_tex_w <<= 1);
	for (GLOBAL_tex_h = 1; GLOBAL_tex_h < GLOBAL_height; GLOBAL_tex_h <<= 1);
 
    // if texture dimensions were changed, realloc the texture memory block (GLOBAL_tex)
	if (GLOBAL_tex_h != oh || GLOBAL_tex_w != ow) {
	    GLOBAL_tex_size = GLOBAL_tex_h * sizeof(rgb_t*) + GLOBAL_tex_h * GLOBAL_tex_w * 3;
		GLOBAL_tex = realloc(GLOBAL_tex, GLOBAL_tex_size);
		// a texture has GLOBAL_tex_h pointers to the begining of each line,
		// followed by the GLOBAL_tex_h lines of the image (bottom to top);
        // each line is sequence of GLOBAL_tex_w*3 bytes (each pixel RGB values)
    }
  
    // update pointers in the beggining of the texture
	for (GLOBAL_tex[0] = (rgb_t *)(GLOBAL_tex + GLOBAL_tex_h), i = 1; i < GLOBAL_tex_h; i++)
		GLOBAL_tex[i] = GLOBAL_tex[i - 1] + GLOBAL_tex_w;
		// uses rgb_t* arithmetic pointer, where each unit corresponds to 3 bytes
}

////////////////////////////////////////////////////////////////////////
void set_texture()
{
	if (GLOBAL_refresh) 
		alloc_tex(); // allocates GLOBAL_tex if necessary and updates GLOBAL_tex_w, GLOBAL_tex_h and GLOBAL_tex_size, 
					 
	calc_mandel(); // calculates the image and stores it in GLOBAL_tex
 
	if (GLOBAL_refresh) { // if the texture was reallocated, create the texture
		glEnable(GL_TEXTURE_2D); // enable 2D textures	
		glBindTexture(GL_TEXTURE_2D, GLOBAL_texture); // bind the texture to the texture object
		glTexImage2D(
			GL_TEXTURE_2D, 0, 3, GLOBAL_tex_w, GLOBAL_tex_h, 
			0, GL_RGB, GL_UNSIGNED_BYTE, GLOBAL_tex[0]); // create the texture
 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // set texture parameters
		render();
	}
}

////////////////////////////////////////////////////////////////////////
void resize(int w, int h)
{
	GLOBAL_width = w;
	GLOBAL_height = h;
 
	glViewport(0, 0, w, h);
	glOrtho(0, w, 0, h, -1, 1);
 
	set_texture();
}

////////////////////////////////////////////////////////////////////////
void mouseclick(int button, int state, int x, int y)
{
	if (state != GLUT_UP) return;
	
	GLOBAL_cx += (x - GLOBAL_width / 2) * GLOBAL_scale;
	GLOBAL_cy -= (y - GLOBAL_height/ 2) * GLOBAL_scale;
 
	switch(button) {
	case GLUT_LEFT_BUTTON: /* zoom in */
		if (GLOBAL_scale > fabs(x) * 1e-16 && GLOBAL_scale > fabs(y) * 1e-16)
			GLOBAL_scale /= 2;
		break;
	case GLUT_RIGHT_BUTTON: /* zoom out */
		GLOBAL_scale *= 2;
		break;
	/* any other button recenters */
	}
	set_texture();
	//print_menu(); // uncomment for convenience; comment for benchmarking
}

////////////////////////////////////////////////////////////////////////
void keypress(unsigned char key, int x, int y)
{
	static int zoomin_x=0, zoomin_y=1; // RUF
	// Ruf: where to start fetching mouse coordinates from GLOBAL_zoomin
	//      (first x coordinate is at GLOBAL_zoomin[0])
	//      (first y coordinate is at GLOBAL_zoomin[1])
	//      (next coordinates are at distance 2: see +=2 at 'z' and 'Z" bellow)
	
	switch(key) {
	case 'q': 
             glFinish();
             glutDestroyWindow(GLOBAL_gwin);
             free(GLOBAL_tex);
             return;
	
	case 27: // Esc
	         GLOBAL_scale = 1./256;
	         GLOBAL_cx = -.6;
	         GLOBAL_cy = 0; 
	         break;
 
	case 'r':
             GLOBAL_color_rotate = (GLOBAL_color_rotate + 1) % 6;
              break;

	case '>':
    case '.':
             GLOBAL_max_iter += 128;
             if (GLOBAL_max_iter > 1 << 15) GLOBAL_max_iter = 1 << 15;
             printf("max iter: %d\n", GLOBAL_max_iter);
             break;
 
	case '<':
    case ',':
             GLOBAL_max_iter -= 128;
             if (GLOBAL_max_iter < 128) GLOBAL_max_iter = 128;
             printf("max iter: %d\n", GLOBAL_max_iter);
             break;

	case 'c':
             GLOBAL_saturation = 1 - GLOBAL_saturation;
             break;
 
	case 's':screen_dump(); return;

	case 'I':GLOBAL_max_iter = 4096; break;
	
	case 'i':GLOBAL_max_iter = 256; break;
	
	case ' ':GLOBAL_invert = !GLOBAL_invert; break;
	
	case 'z':// simulate one mouse click in order to dive one time in zoomin
             GLOBAL_refresh=1;
             mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]);
             zoomin_x+=2; zoomin_y+=2;
             break;

	case 'Z':// simulate many mouse clicks in order to dive fully in zoomin
             GLOBAL_refresh=1; // use 0 to avoid refreshing all but the last one
             for (zoomin_x=0, zoomin_y=1; zoomin_x < GLOBAL_zoomin_num_pairs; zoomin_x+=2, zoomin_y +=2) {
                 if (zoomin_x == GLOBAL_zoomin_num_pairs-2) GLOBAL_refresh=1;
                 mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]); 					
             }

             // simulate case 's'
             keypress('s', -1, -1);
              
             // simulate case 'q'
             keypress('q', -1, -1);
             return;
	}
	set_texture();
	print_menu();
}

////////////////////////////////////////////////////////////////////////
void render()
{
	double	x = (double)GLOBAL_width /GLOBAL_tex_w,
		    y = (double)GLOBAL_height/GLOBAL_tex_h;
 
	glClear(GL_COLOR_BUFFER_BIT);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
	glBindTexture(GL_TEXTURE_2D, GLOBAL_texture);
 
	glBegin(GL_QUADS);
 
	glTexCoord2f(0, 0); glVertex2i(0, 0);
	glTexCoord2f(x, 0); glVertex2i(GLOBAL_width, 0);
	glTexCoord2f(x, y); glVertex2i(GLOBAL_width, GLOBAL_height);
	glTexCoord2f(0, y); glVertex2i(0, GLOBAL_height);
 
	glEnd();
 
	glFlush();
	glFinish();
}

////////////////////////////////////////////////////////////////////////
// initializes the graphics
// complexity: O(1)
// note: this function is called only once
// 	 (at the beginning of the program)

void init_gfx(int *c, char **v)
{
	glutInit(c, v); // initialize GLUT library state c and v are command line arguments
	glutInitDisplayMode(GLUT_RGB); // set display mode to RGB color mode and single buffering 
	glutInitWindowSize(GLOBAL_window_width, GLOBAL_window_height); // set initial window size to 1024x768 pixels 
 
	GLOBAL_gwin = glutCreateWindow("Mandelbrot"); // create window with title "Mandelbrot" and return window id
	glutDisplayFunc(render); // set display callback for current window to render function, which is called when window needs to be redrawn 
 
	glutKeyboardFunc(keypress); // set keyboard callback for current window to keypress function, which is called when a key is pressed
	glutMouseFunc(mouseclick); // set mouse callback for current window to mouseclick function, which is called when a mouse button is pressed
	glutReshapeFunc(resize);  // set reshape callback for current window to resize function, which is called when window is resized
	glGenTextures(1, &GLOBAL_texture); // generate a texture name and store it in GLOBAL_texture
	set_texture(); // set texture parameters
}

////////////////////////////////////////////////////////////////////////
// main function 
// responsible for initializing the graphics and starting the main loop
// complexity: O(1)

int main(int c, char **v){
	init_gfx(&c, v);
	print_menu();
	glutMainLoop();	
	return 0;
}
