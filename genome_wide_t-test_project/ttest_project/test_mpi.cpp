#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <numeric>      // std::iota
#include <algorithm>
#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>
#include <tuple>
#include <utility>
using namespace std;

//#define SIZE 4532
//#define MAX_RAND 100
//float D[SIZE][8], C[SIZE][52], T[SIZE][1];
//float **D; float **C; float **T;

float **alloc_2d_float(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}
// Function to find mean.
float Mean(float arr[], int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + arr[i];
    return sum / n;
}
 
// Function to find standard
// deviation of given array.
float standardDeviation(float arr[], int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + (arr[i] - Mean(arr, n)) *
                    (arr[i] - Mean(arr, n));
 
    return sqrt(sum / (n - 1));
}
 
// Function to find t-test of
// two set of statistical data.
float tTest(float arr1[], int n,
            float arr2[], int m)
{
    float mean1 = Mean(arr1, n);
    float mean2 = Mean(arr2, m);
    float sd1 = standardDeviation(arr1, n);
    float sd2 = standardDeviation(arr2, m);
 
    // Formula to find t-test
    // of two set of data.
    float t_test = (mean1 - mean2) / sqrt((sd1 * sd1)
                              / n + (sd2 * sd2) / m);
    return t_test;
}

#define SIZE 4000       /* Size of matrices */
#define MAX_RAND 10
#define colSIZE 52

float A[SIZE][colSIZE], B[SIZE][8];//D[SIZE][SIZE],E[SIZE][1];
float C[SIZE][1];
int k = 5;

float *returnKItems(float *arr, int k) {
    float *arrK = new float[k];
    int nRows = sizeof (arr) / sizeof (arr[0]);
    random_device rd;
    mt19937 g(rd());
    // Use a different seed value so that we don't get
    // same result each time we run this program
    srand(time(0));
    int j = 0;
    do {
        for (int i = 0; i < k; i++) {
        arrK[i] = arr[k];
        }
        j++;
    } while ((next_permutation(arr, arr + nRows)) & (j < 1));
    //shuffle(arr, arr+(), g);
    
    return arrK;
}

float *replaceItems(float *arr, float *oldArr, float *newArr)
{
    for (int i=0; i < sizeof (oldArr) / sizeof (oldArr[0]); i++)
    {
        replace(arr, arr + (sizeof (arr) / sizeof (arr[0])), oldArr[i], newArr[i]);
    }
    return arr;
}

float *returnRandomPerm(float *arr1, float *arr2, int k)
{
    float *randomTtests = new float[1000];
    int arr1Size = sizeof (arr1) / sizeof (arr1[0]);
    int arr2Size = sizeof (arr2) / sizeof (arr2[0]);
    for (int m = 0; m < 1000; m++) {
        float *randArr = returnKItems(arr1, k); //select k items - 5
        float *randArr2 = returnKItems(arr2, k);
        float *newArr = replaceItems(arr1, randArr, randArr2);
        float *newArr2 = replaceItems(arr2, randArr2, randArr);
        randomTtests[m] = tTest(newArr, arr1Size, newArr2, arr2Size);
    }
    return randomTtests;
}

void fill_matrix(float m[SIZE][SIZE])
{
    //static int n=0;
    int i, j;

    printf("\n*****************************\n");
    for (i=0; i<SIZE; i++)
    {
        for (j=0; j<SIZE; j++){     
            m[i][j] = rand() % MAX_RAND;
            printf("%2f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n*****************************\n");
}

void fill_vector(float m[SIZE][1])
{
    //static int n=0;
    int i, j;

    printf("\n*****************************\n");
    for (i=0; i<SIZE; i++)
    {
        for (j=0; j<1; j++){     
            m[i][j] = rand() % MAX_RAND;
            printf("%2f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n*****************************\n");
}


void print_matrix(float m[SIZE][SIZE])
{
    int i, j = 0;
    for (i=0; i<SIZE; i++) {
        printf("\n\t| ");
        for (j=0; j<SIZE; j++)
            printf("%2f ", m[i][j]);
        printf("|");
    }
}

void print_vector(float m[SIZE][1])
{
    int i, j = 0;
    for (i=0; i<SIZE; i++) {
        printf("\n\t| ");
        for (j=0; j<1; j++)
            printf("%2f ", m[i][j]);
        printf("|");
    }
}


int main(int argc, char *argv[])
{
    int myrank, P, from, to, i, j, k;
    /*float **A = alloc_2d_float(SIZE, colSIZE);
    float **B = alloc_2d_float(SIZE, 8);
    float **C = alloc_2d_float(SIZE, 1);*/
    //  int tag = 666;      /* any value will do */
    //  MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   /* who am i */
    MPI_Comm_size(MPI_COMM_WORLD, &P); /* number of processors */


    if (SIZE%P!=0) {
        if (myrank==0) printf("Matrix size not divisible by number of processors\n");
        MPI_Finalize();
        exit(-1);
    }

    from = myrank * SIZE/P;
    to = ((myrank+1) * SIZE/P);

 /* Process 0 fills the input matrices and broadcasts them to the rest */
    /* (actually, only the relevant stripe of A is sent to each process) */

    if (myrank==0) {


        //static int n=0;
        int i, j;

        printf("\n*********FILL A Matrix: 64*64 ********************\n");
        for (i=0; i<SIZE; i++)
        {
            for (j=0; j<colSIZE; j++){     
                A[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                //rand() % MAX_RAND;        
                printf("%f ", A[i][j]);
            }
            printf("\n");

            printf("\n*****************************\n");
        }
        printf("\n**********Fill B Matrix 64 * 8 *******************\n");
        for (i=0; i<SIZE; i++)
        {
            for (j=0; j<8; j++){     
                B[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                //rand() % MAX_RAND;        
                printf("%f ", B[i][j]);
            }
            //C[i][0]=0;
            printf("\n");

            printf("\n*****************************\n");
        }
        //fill_vector(B);
    }

    //int s=SIZE*SIZE/P;
    // printf("computing slice %d (from row %d to %d)\n", myrank, from, to-1);
    MPI_Bcast (B, SIZE*8, MPI_INT, 0, MPI_COMM_WORLD);
    //    printf("\n\n%d",s);
    //print_vector(s);
    //printf("\n\n");
    if(myrank==0){
        MPI_Scatter (&A[0][0], SIZE*colSIZE/P, MPI_INT, MPI_IN_PLACE, SIZE*colSIZE/P, MPI_INT, 0, MPI_COMM_WORLD);
    }else{
        MPI_Scatter (&A[0][0], SIZE*colSIZE/P, MPI_INT, &A[from][0], SIZE*colSIZE/P, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("computing slice %d (from row %d to %d)\n", myrank, from, to-1);
    for (i=from; i<to; i++) {
        C[i][0]=0;
        float subArr[colSIZE];
        for (k=0; k<colSIZE; k++){
            subArr[k] = A[i][k];
        }
        float subArr2[8];
        for (int v=0; v<8; v++){
            subArr2[v] = B[i][v];
        }
        float tTst= tTest(subArr, colSIZE, subArr2, 8);
        //C[i][0] = tTst;
        float randomPerms[1000];
        for (int n = 0; n < 1000; n++) {
            float *randArr = returnKItems(subArr, k); //select k items - 5
            float *randArr2 = returnKItems(subArr2, k);
            float *newArr = replaceItems(subArr, randArr, randArr2);
            float *newArr2 = replaceItems(subArr2, randArr2, randArr);
            randomPerms[n] = tTest(newArr, colSIZE, newArr2, 8);
        }
        //float *randomPerms = returnRandomPerm(subArr, subArr2, k);
        //int randIndex = rand() % 1000;
        float dScore = Mean(randomPerms, 1000);
        //sizeof (randomPerms) / sizeof (randomPerms[0]);
        //(Mean(randomTtestVec, 1000)/standardDeviation(randomTtestVec, 1000));
        C[i][0] = dScore;
        fill(randomPerms, randomPerms+1000, 0);
        fill(subArr, subArr+colSIZE, 0);
        fill(subArr2, subArr2+8, 0);
        
    }
    
    if(myrank==0){
         MPI_Gather (MPI_IN_PLACE, SIZE/P, MPI_INT, &C[0][0], SIZE/P, MPI_INT, 0, MPI_COMM_WORLD);
    }else{
         MPI_Gather (&C[from][0], SIZE/P, MPI_INT, &C[0][0], SIZE/P, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (myrank==0) {
        printf("\n\n");
        {
            int i, j = 0;
            for (i=0; i<SIZE; i++) {
                printf("Gene Index:  %d ", i);
                printf("\nD score: %f ", C[i][0]);
                printf("\n");
            }
        }

        printf("\n\n");
        //   print_matrix(D);
        //printf("\n\n\t       * \n");
        //print_vector(B);
        //printf("\n\n\t       = \n");
        //print_vector(C);
        //printf("\n\n");
        //   print_vector(E);
        //   printf("\n\n");
    }
    MPI_Finalize();
    return 0;
}
