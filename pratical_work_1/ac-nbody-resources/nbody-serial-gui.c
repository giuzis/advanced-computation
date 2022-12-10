#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <graphics.h>

#define CONST_GravConstant 0.01
 
typedef struct{
	double x,y,z;
}vector;
 
int GLOBAL_numBodies=15; // default number of bodies
int GLOBAL_numSteps=30000; // default number of time steps (1 million)

int GLOBAL_windowWidth=800; // default window width is 800 pixels
int GLOBAL_windowHeight=800; // default window height is 800 pixels


const char GLOBAL_useGUI=1; // use GUI by default
int GLOBAL_pixel_or_circle=0; // 1 = bodies will be pixels; 0 = bodies will be circles
int GLOBAL_circle_radius_add=3; // circle radius is equal to the body mass plus this value
char GLOBAL_clearDevice=1; // pause at the end; default is yes (1)
char GLOBAL_finalPause=1; // pause at the end; default is yes (1)


double *GLOBAL_masses;
vector *GLOBAL_positions;
vector *GLOBAL_velocities;
vector *GLOBAL_accelerations;

//////////////////////////////////////////////////////////////////////// 
vector addVectors(vector a,vector b){
	vector c = {a.x+b.x,a.y+b.y,a.z+b.z};
 
	return c;
}

//////////////////////////////////////////////////////////////////////// 
vector scaleVector(double b,vector a){
	vector c = {b*a.x,b*a.y,b*a.z};
 
	return c;
}
 
////////////////////////////////////////////////////////////////////////
vector subtractVectors(vector a,vector b){
	vector c = {a.x-b.x,a.y-b.y,a.z-b.z};
 
	return c;
}

////////////////////////////////////////////////////////////////////////
double mod(vector a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

////////////////////////////////////////////////////////////////////////
void initSystemFromRandom(){
	int i;		

	GLOBAL_masses = (double*)malloc(GLOBAL_numBodies*sizeof(double));
	GLOBAL_positions = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
	GLOBAL_velocities = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
	GLOBAL_accelerations = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
 
	for(i=0;i<GLOBAL_numBodies;i++){
		GLOBAL_masses[i]=rand()%GLOBAL_numBodies;
		GLOBAL_positions[i].x=rand()%GLOBAL_windowWidth;
		GLOBAL_positions[i].y=rand()%GLOBAL_windowHeight;
		GLOBAL_positions[i].z=0;
		GLOBAL_velocities[i].x=GLOBAL_velocities[i].y=GLOBAL_velocities[i].z=0;
	}
}

////////////////////////////////////////////////////////////////////////
void showSystem(){
	int i;

	for(i=0;i<GLOBAL_numBodies;i++) {
		fprintf(stderr,"Body %d : %lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i,\
		GLOBAL_positions[i].x,GLOBAL_positions[i].y,GLOBAL_positions[i].z,\
		GLOBAL_velocities[i].x,GLOBAL_velocities[i].y,GLOBAL_velocities[i].z);
	}
}

////////////////////////////////////////////////////////////////////////
void validateSystem(){
	int i;
	
	for(i=0;i<GLOBAL_numBodies;i++) {
		if (isnan(GLOBAL_positions[i].x) || isnan(GLOBAL_positions[i].y) || isnan(GLOBAL_positions[i].z) || \
		    isnan(GLOBAL_velocities[i].x) || isnan(GLOBAL_velocities[i].y) || isnan(GLOBAL_velocities[i].z) ) {
			fprintf(stderr,"NAN Body %d : %lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i,\
			GLOBAL_positions[i].x,GLOBAL_positions[i].y,GLOBAL_positions[i].z,\
			GLOBAL_velocities[i].x,GLOBAL_velocities[i].y,GLOBAL_velocities[i].z);
			exit(1);
		}
	}
}

//////////////////////////////////////////////////////////////////////// 
void computeAccelerations(){
	int i,j;
 
	for(i=0;i<GLOBAL_numBodies;i++){
		GLOBAL_accelerations[i].x = 0;
		GLOBAL_accelerations[i].y = 0;
		GLOBAL_accelerations[i].z = 0;
		for(j=0;j<GLOBAL_numBodies;j++){
			if(i!=j){
				GLOBAL_accelerations[i] = addVectors(GLOBAL_accelerations[i],scaleVector(CONST_GravConstant*GLOBAL_masses[j]/pow(mod(subtractVectors(GLOBAL_positions[i],GLOBAL_positions[j])),3),subtractVectors(GLOBAL_positions[j],GLOBAL_positions[i])));
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////// 
void computeVelocities(){
	int i;

	for(i=0;i<GLOBAL_numBodies;i++)
		GLOBAL_velocities[i] = addVectors(GLOBAL_velocities[i],GLOBAL_accelerations[i]);
}

//////////////////////////////////////////////////////////////////////// 
void computePositions(){
	int i;

	for(i=0;i<GLOBAL_numBodies;i++)
		GLOBAL_positions[i] = addVectors(GLOBAL_positions[i],addVectors(GLOBAL_velocities[i],scaleVector(0.5,GLOBAL_accelerations[i])));
}

////////////////////////////////////////////////////////////////////////
void resolveCollisions(){
	int i,j;

	for(i=0;i<GLOBAL_numBodies-1;i++)
		for(j=i+1;j<GLOBAL_numBodies;j++){
			if(GLOBAL_positions[i].x==GLOBAL_positions[j].x && GLOBAL_positions[i].y==GLOBAL_positions[j].y && GLOBAL_positions[i].z==GLOBAL_positions[j].z){
				vector temp = GLOBAL_velocities[i];
				GLOBAL_velocities[i] = GLOBAL_velocities[j];
				GLOBAL_velocities[j] = temp;
			}
		}
}

//////////////////////////////////////////////////////////////////////// 
void simulate(){
	computeAccelerations();
	computePositions();
	computeVelocities();
	resolveCollisions();
	validateSystem();
}

//////////////////////////////////////////////////////////////////////// 
void plotOrbits(){
	int i;
 
	for(i=0;i<GLOBAL_numBodies;i++) {		
		if (GLOBAL_pixel_or_circle)
			putpixel(GLOBAL_positions[i].x,GLOBAL_positions[i].y, i%15 + 1);
		else {
			setcolor(i%15 + 1);
			circle(GLOBAL_positions[i].x, GLOBAL_positions[i].y, GLOBAL_masses[i]+GLOBAL_circle_radius_add);
		}
    }
}

////////////////////////////////////////////////////////////////////////
int main(int argC,char* argV[])
{
	int i; char s[128];		
	const unsigned int randSeed=10; // default value; do not change
				
	if (argC>=3) {		
		GLOBAL_numBodies=atoi(argV[1]);
		GLOBAL_numSteps=atoi(argV[2]);
		
		if(argC>=5) {
			GLOBAL_windowWidth=atoi(argV[3]);
			GLOBAL_windowHeight=atoi(argV[4]);
		}
	}
	else if (argC!=1) {
		printf("Usage : %s [numBodies numSteps [windowWidth windowHeight]]\n",argV[0]);
		exit(1);
	}
	
	if (GLOBAL_useGUI)
		initwindow(GLOBAL_windowWidth,GLOBAL_windowHeight);
	
	srand(randSeed);
	initSystemFromRandom();
	//showSystem();

	for(i=0;i<GLOBAL_numSteps;i++){
		
		if (GLOBAL_useGUI) {
			if (GLOBAL_clearDevice)
				cleardevice(); // remove the objects track
			sprintf (s, "%d", i+1); 
			setcolor (WHITE); outtextxy (0, 0, s);			
		}
		
		//printf("%d\n", i);			
		
		simulate();	
		//showSystem();
		
		if (GLOBAL_useGUI) {
			plotOrbits();												
			refresh();					
			if (kbhit()) break;
		}

	}

	// dump final result
	showSystem();
	
	if (GLOBAL_finalPause) {
		sprintf (s, "simulation stoped or finished; press ENTER to quit"); 
		if (GLOBAL_useGUI) {			
			setcolor (WHITE); outtextxy (0, 10, s); refresh ();					
			while (!kbhit());
			closegraph ();
		 }
		 else
			getchar();			
	}
		
	return 0;
}
 
