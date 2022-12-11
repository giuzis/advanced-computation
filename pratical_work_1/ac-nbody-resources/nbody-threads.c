#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define CONST_GravConstant 0.01
#define NUM_THREADS 8
 
typedef struct{
	double x, y, z;
}vector;
 
int GLOBAL_numBodies=15; // default number of bodies
int GLOBAL_numSteps=30000; // default number of time steps (1 million)
int TOTAL_WORK = GLOBAL_numBodies*(GLOBAL_numBodies - 1)/(2); 
int WORK_PER_THREAD = TOTAL_WORK/(N_THREADS); 

int GLOBAL_windowWidth=800; // default window width is 800 pixels
int GLOBAL_windowHeight=800; // default window height is 800 pixels

double *GLOBAL_masses;
vector *GLOBAL_positions;
vector *GLOBAL_velocities;
vector *GLOBAL_accelerations;
vector **GLOBAL_partial_accelerations;

//////////////////////////////////////////////////////////////////////// 
typedef struct{
	int start;
	int end;
}range;

pthread_t threads[NUM_THREADS];
range slices[NUM_THREADS];

// create barrier to global aceleration
pthread_barrier_t GLOBAL_barrier;

int GLOBAL_WORK_PER_THREAD = 0;

void updateWorkPerThread(void){
	int work_amount = ((GLOBAL_numBodies*(GLOBAL_numBodies - 1))/2);
	GLOBAL_WORK_PER_THREAD = work_amount/(NUM_THREADS);
	printf("Work amount: %d, Work per thread: %d\n", work_amount, GLOBAL_WORK_PER_THREAD);
}

int sliceSize(int start){
	int b = 2*(GLOBAL_numBodies - start) - 1;
    int range = (b - sqrt(b * b - 8 * GLOBAL_WORK_PER_THREAD )) / 2;
	int result = start + range + (range == 0);
	if ((result > GLOBAL_numBodies - 1) || (range < 0))
		return GLOBAL_numBodies - 1;
	return result;
}

void calculateSlices(void){
	for(int i = 0, start = 0; i < NUM_THREADS; i++){
		slices[i].start = start;
		slices[i].end = (start = sliceSize(start));
	}
	slices[NUM_THREADS - 1].end = GLOBAL_numBodies - 1;	
	for(int i = 0; i < NUM_THREADS; i++)
		printf("Slice %d: %d - %d\n", i, slices[i].start, slices[i].end);
}
//////////////////////////////////////////////////////////////////////// 
vector addVectors(vector a, vector b){
	return (vector){a.x + b.x, a.y + b.y, a.z + b.z};
}
//////////////////////////////////////////////////////////////////////// 
vector invertVector(vector a){
	return (vector){-a.x, -a.y, -a.z};
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
	GLOBAL_partial_accelerations = (vector**)malloc(GLOBAL_numBodies*sizeof(vector*));
	
	for(int i = 0; i < GLOBAL_numBodies; i++){
		GLOBAL_masses[i] = rand() % GLOBAL_numBodies;
		GLOBAL_positions[i] = (vector){rand() % GLOBAL_windowWidth, rand() % GLOBAL_windowHeight, 0};
		GLOBAL_velocities[i] = (vector){0, 0, 0};
	}
}
////////////////////////////////////////////////////////////////////////
void initSystemMatrix(int start, int end){
	for(int i = start; i < end; i++)
		GLOBAL_partial_accelerations[i] = (vector*)malloc(GLOBAL_numBodies*sizeof(vector));
}
////////////////////////////////////////////////////////////////////////
void showSystem(){
	for(int i = 0; i < GLOBAL_numBodies; i++) {
		printf("Body %d : %lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i,
		GLOBAL_positions[i].x, GLOBAL_positions[i].y, GLOBAL_positions[i].z,
		GLOBAL_velocities[i].x, GLOBAL_velocities[i].y, GLOBAL_velocities[i].z);
	}
}

////////////////////////////////////////////////////////////////////////
void validateSystem(int start, int end){
	for(int i = start; i < end; i++) {
		if (isnan(GLOBAL_positions[i].x) || isnan(GLOBAL_positions[i].y) || isnan(GLOBAL_positions[i].z) || 
		    isnan(GLOBAL_velocities[i].x) || isnan(GLOBAL_velocities[i].y) || isnan(GLOBAL_velocities[i].z) ) {
			printf("NAN Body %d : %lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i,\
			GLOBAL_positions[i].x,GLOBAL_positions[i].y,GLOBAL_positions[i].z,\
			GLOBAL_velocities[i].x,GLOBAL_velocities[i].y,GLOBAL_velocities[i].z);
			exit(1);
		}
	}
}
//////////////////////////////////////////////////////////////////////// 
void computePartialAccelerations(int start, int end){
	for(int i = start; i < end; i++)
		for(int j = i + 1; j < GLOBAL_numBodies; j++){
			vector distance = subtractVectors(GLOBAL_positions[j], GLOBAL_positions[i]);
			double distance_mod3 = pow(mod(distance), 3);
			
			vector acceleration_i = scaleVector(( CONST_GravConstant * GLOBAL_masses[j]) / distance_mod3, distance);
			vector acceleration_j = scaleVector((-CONST_GravConstant * GLOBAL_masses[i]) / distance_mod3, distance);
			
			GLOBAL_partial_accelerations[i][j] = acceleration_i;
			GLOBAL_partial_accelerations[j][i] = acceleration_j;
		}
}
//////////////////////////////////////////////////////////////////////// 
void computeAccelerations(int start, int end){
	for(int i = start; i < end; i++){
		GLOBAL_accelerations[i] = (vector){0, 0, 0};
		for(int j = 0; j < GLOBAL_numBodies; j++)
			if(i != j)
				GLOBAL_accelerations[i] = addVectors(GLOBAL_accelerations[i], GLOBAL_partial_accelerations[i][j]);			
	}

	for(int j = 0; j < GLOBAL_numBodies; j++){
		GLOBAL_accelerations[j] = addVectors(GLOBAL_accelerations[j], local_accelerations[j]);
	}

}
//////////////////////////////////////////////////////////////////////// 
void computeVelocities(int start, int end){
	for(int i = start; i < end; i++)
		GLOBAL_velocities[i] = addVectors(GLOBAL_velocities[i], GLOBAL_accelerations[i]);
}

//////////////////////////////////////////////////////////////////////// 
void computePositions(int start, int end){
	for(int i = start; i < end; i++)
		GLOBAL_positions[i] = addVectors(GLOBAL_positions[i],
			addVectors(GLOBAL_velocities[i], scaleVector(0.5, GLOBAL_accelerations[i])));
}

////////////////////////////////////////////////////////////////////////
void resolveCollisions(int start, int end){
	for(int i = start; i < end; i++)
		for(int j = i + 1; j < GLOBAL_numBodies; j++){
			if(GLOBAL_positions[i].x == GLOBAL_positions[j].x && GLOBAL_positions[i].y == GLOBAL_positions[j].y && GLOBAL_positions[i].z ==GLOBAL_positions[j].z){
				vector temp = GLOBAL_velocities[i];
				GLOBAL_velocities[i] = GLOBAL_velocities[j];
				GLOBAL_velocities[j] = temp;
			}
		}
}

//////////////////////////////////////////////////////////////////////// 
void* simulate(void* arg){
	long tid = (long)arg;
	
	int start_range = slices[tid].start; 
	int end_range = slices[tid].end;

	int work = GLOBAL_numBodies/NUM_THREADS, start_normal = -1, end_normal = -1;
	if (work){
		start_normal = tid*(work);
		end_normal = tid == NUM_THREADS - 1? GLOBAL_numBodies : (start_normal + work);
	}
	else if(tid < GLOBAL_numBodies){
		start_normal = tid;
		end_normal = tid + 1;
	}	
	printf("Thread %ld : %d - %d\n", tid, start_normal, end_normal);
	initSystemMatrix(start_normal, end_normal);
	pthread_barrier_wait(&GLOBAL_barrier);
	for(int i = 0; i < GLOBAL_numSteps; i++){
		computePartialAccelerations(start_range, end_range);

		pthread_barrier_wait(&GLOBAL_barrier);
		computeAccelerations(start_normal, end_normal);
		computePositions(start_normal, end_normal);
		computeVelocities(start_normal, end_normal);

		pthread_barrier_wait(&GLOBAL_barrier);
		resolveCollisions(start_range, end_range);
		pthread_barrier_wait(&GLOBAL_barrier);
		validateSystem(start_normal, end_normal);
		pthread_barrier_wait(&GLOBAL_barrier);
	}
	pthread_exit(NULL);
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

	TOTAL_WORK = GLOBAL_numBodies*(GLOBAL_numBodies - 1)/(2); 
	WORK_PER_THREAD = TOTAL_WORK/(N_THREADS); 
		
	srand(randSeed);
	initSystemFromRandom();
	//showSystem();
	updateWorkPerThread();
	calculateSlices();
	
	int ret = pthread_barrier_init(&GLOBAL_barrier, NULL, NUM_THREADS);
    if (ret){
        printf("ERROR; return code from pthread_barrier_init() is %d\n", ret);
        exit(ret);
    }
	
	for(long tid = 0; tid < NUM_THREADS; tid++){
		int ret = pthread_create(&threads[tid], NULL, simulate, (void*)tid);
		if(ret){
			printf("ERROR; return code from pthread_create() is %d\n", ret);
			exit(ret);
		}
	}
	for(long tid = 0; tid < NUM_THREADS; tid++){
		int ret = pthread_join(threads[tid], NULL);
		if(ret){
			printf("ERROR; return code from pthread_join() is %d\n", ret);
			exit(ret);
		}
	}				
	showSystem();
	return 0;
}
 

