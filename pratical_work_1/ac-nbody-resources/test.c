#include<stdio.h>

typedef struct{
	double x,y,z;
}vector;

vector addVectors(vector a, vector b){
	return (vector){a.x + b.x, a.y + b.y, a.z + b.z};
}

vector* addVectorsP(vector *a, vector b){
	*a = (vector){a->x + b.x, a->y + b.y, a->z + b.z};
    return a;
}

vector invertVector(vector a){
    return (vector){-a.x, -a.y, -a.z};
}

vector scaleVector(double b, vector a){
    return (vector){b * a.x, b * a.y, b * a.z};
}

vector subtractVectors(vector a,vector b){
    return (vector){a.x - b.x, a.y - b.y, a.z - b.z};
}

double mod(vector a){
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double mod3(vector a){
    double mod2 = (a.x * a.x + a.y * a.y + a.z * a.z); 
    return sqrt(mod2) * mod2;
}


int main(){
    vector a = {1,2,3};
    for (int i = 0; i < 1000000000; i++)
    // a = addVectors(a, (vector){1,2,3});
    addVectorsP(&a, (vector){1,2,3});
    printf("%lf %lf %lf\n", a.x, a.y, a.z);
    return 0;
}