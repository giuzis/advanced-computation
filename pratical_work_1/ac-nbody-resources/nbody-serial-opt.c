#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CONST_GravConstant 0.01
 
typedef struct{
	double x, y, z;
}vector;
 
int GLOBAL_numBodies=15; // default number of bodies
int GLOBAL_numSteps=30000; // default number of time steps (1 million)

int GLOBAL_windowWidth=800; // default window width is 800 pixels
int GLOBAL_windowHeight=800; // default window height is 800 pixels

double *GLOBAL_masses;
vector *GLOBAL_positions;
vector *GLOBAL_velocities;
vector *GLOBAL_accelerations;

//////////////////////////////////////////////////////////////////////// 
vector addVectors(vector a, vector b){
	return (vector){a.x + b.x, a.y + b.y, a.z + b.z};
}
//////////////////////////////////////////////////////////////////////// 
vector scaleVector(double b, vector a){
	return (vector){b * a.x, b * a.y, b * a.z};
}
////////////////////////////////////////////////////////////////////////
vector subtractVectors(vector a,vector b){
	return (vector){a.x - b.x, a.y - b.y, a.z - b.z};
}
////////////////////////////////////////////////////////////////////////
double mod(vector a){
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
////////////////////////////////////////////////////////////////////////
void initSystemFromRandom(){
	GLOBAL_masses = (double*)malloc(GLOBAL_numBodies*sizeof(double));
	GLOBAL_positions = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
	GLOBAL_velocities = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
	GLOBAL_accelerations = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
	
	for(int i = 0; i < GLOBAL_numBodies; i++){
		GLOBAL_masses[i] = rand() % GLOBAL_numBodies;
		GLOBAL_positions[i] = (vector){rand() % GLOBAL_windowWidth, rand() % GLOBAL_windowHeight, 0};
		GLOBAL_velocities[i] = (vector){0, 0, 0};
	}
}

////////////////////////////////////////////////////////////////////////
void showSystem(){
	for(int i = 0; i < GLOBAL_numBodies; i++) {
		fprintf(stderr,"Body %d : %lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i,
		GLOBAL_positions[i].x, GLOBAL_positions[i].y, GLOBAL_positions[i].z,
		GLOBAL_velocities[i].x, GLOBAL_velocities[i].y, GLOBAL_velocities[i].z);
	}
}

////////////////////////////////////////////////////////////////////////
void validateSystem(){
	for(int i = 0; i < GLOBAL_numBodies; i++) {
		if (isnan(GLOBAL_positions[i].x) || isnan(GLOBAL_positions[i].y) || isnan(GLOBAL_positions[i].z) || 
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
	for(int i = 0; i < GLOBAL_numBodies; i++)
		GLOBAL_accelerations[i] = (vector){0, 0, 0};
	
	for(int i = 0; i < GLOBAL_numBodies - 1; i++){
		for(int j = i + 1; j < GLOBAL_numBodies; j++){
			vector distance = subtractVectors(GLOBAL_positions[j], GLOBAL_positions[i]);
			double distance_mod3 = pow(mod(distance), 3);

			GLOBAL_accelerations[i] = addVectors(GLOBAL_accelerations[i], scaleVector(CONST_GravConstant * GLOBAL_masses[j] / distance_mod3, distance));
			GLOBAL_accelerations[j] = addVectors(GLOBAL_accelerations[j], scaleVector(CONST_GravConstant * -GLOBAL_masses[i] / distance_mod3, distance));
		}
	}
}


//////////////////////////////////////////////////////////////////////// 
void computeVelocities(){
	for(int i = 0; i < GLOBAL_numBodies; i++)
		GLOBAL_velocities[i] = addVectors(GLOBAL_velocities[i], GLOBAL_accelerations[i]);
}

//////////////////////////////////////////////////////////////////////// 
void computePositions(){
	for(int i = 0; i < GLOBAL_numBodies; i++)
		GLOBAL_positions[i] = addVectors(GLOBAL_positions[i], addVectors(GLOBAL_velocities[i], scaleVector(0.5, GLOBAL_accelerations[i])));
}

////////////////////////////////////////////////////////////////////////
void resolveCollisions(){
	for(int i = 0; i < GLOBAL_numBodies - 1; i++)
		for(int j = i + 1; j < GLOBAL_numBodies; j++){
			if(GLOBAL_positions[i].x == GLOBAL_positions[j].x && GLOBAL_positions[i].y == GLOBAL_positions[j].y && GLOBAL_positions[i].z ==GLOBAL_positions[j].z){
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
int main(int argC,char* argV[]){
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
		
	srand(randSeed);
	initSystemFromRandom();
	//showSystem();

	for(int i = 0; i < GLOBAL_numSteps; i++){
		//printf("%d\n", i);			
		simulate();	
		//showSystem();		
	}
	// dump final result
	showSystem();
		
	return 0;
}
 
