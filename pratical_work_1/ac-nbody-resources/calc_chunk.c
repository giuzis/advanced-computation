#include <stdio.h>
#include <math.h>

#define N_BODIES 4000
#define N_THREADS 8
#define WORK_AMOUNT (N_BODIES*(N_BODIES - 1)/(2))
#define CHUNK_SIZE (WORK_AMOUNT/( N_THREADS))

int get_X(int i, int ntread){
    // int a = -1;
    // int b = 2*N_BODIES - 1 - 2*i;
    // int c = -2*CHUNK_SIZE;
    // int x = (-b + sqrt(b * b - 4 * a * c ))/(2 * a);
    int b = 2*N_BODIES - 1 - 2*i;
    int x = (b - sqrt(b * b - 8 * CHUNK_SIZE )) / 2 ;
    
    printf("-----------------------------------------\nX: %d\n", x);
    
    return ntread == N_THREADS - 1 ? N_BODIES - 1 : x + i ;
}


int tread(int i0, int i1, int tid){
    int sum = 0;
    for(int i = i0; i < i1; i++)
        for(int j = i + 1; j < N_BODIES; j++)
            sum++;
    printf("Thread %d (%d  %d), sum: %d\n", tid, i0, i1, sum);
    return sum;
}

int main(){
    printf("WORK_AMOUNT: %d\n", WORK_AMOUNT);
    printf("CHUNK_SIZE: %d\n", CHUNK_SIZE);
    int i = 0, i1, total = 0, tid = 0;
    do{
        i1 = get_X(i, tid);
        total += tread(i, i1, tid++);
        i = i1;
    }while(i1 < N_BODIES - 2);
    printf("-----------------------------------------\nTotal: %d\n", total);
    return 0;
}
