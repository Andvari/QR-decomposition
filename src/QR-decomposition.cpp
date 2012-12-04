//============================================================================
// Name        : QR-decomposition.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdio.h>
#include <malloc.h>
#include <math.h>

using namespace std;

void  print_matrix(char *, float **, int, int);
void  print_vector(char *, float *, int);
void  makeHouseholder(float **, float **, int);
float makeNorma(float **, int, int);
float sign(float);
void mult_matrix(float **, float **, float **, int);
void makeSolve(float **, float *, float *, int);

int main(void) {
	FILE *f;
	float **A;
	float **Q;
	float **Qt;
	float **R;
	float *b;
	float *y;
	float *x;
	int i, j;
	int N;

	f = fopen("matrix.txt", "r");
	fscanf(f, "%d", &N);

	A = (float **)malloc(N*sizeof(float *));
	Q = (float **)malloc(N*sizeof(float *));
	Qt = (float **)malloc(N*sizeof(float *));
	R = (float **)malloc(N*sizeof(float *));
	b = (float *)malloc(N*sizeof(float *));
	y = (float *)malloc(N*sizeof(float *));
	x = (float *)malloc(N*sizeof(float *));

	for(i=0 ; i<N; i++){
		A[i] = (float *)malloc(N * sizeof(float));
		Q[i] = (float *)malloc(N * sizeof(float));
		Qt[i] = (float *)malloc(N * sizeof(float));
		R[i] = (float *)malloc(N * sizeof(float));
	}

	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			fscanf(f, "%f", &A[i][j]);
		}
	}

	for(j=0; j<N; j++){
		fscanf(f, "%f", &b[j]);
	}

	print_matrix("A", A, N, N);

	print_vector("b", b, N);

	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			R[i][j] = A[i][j];
		}
	}

	makeHouseholder(R, Q, N);

	print_matrix("R", R, N, N);

	print_matrix("Q", Q, N, N);

	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			Qt[i][j] = Q[j][i];
		}
	}

	for(i=0; i<N; i++){
		y[i] = 0;
		for(j=0; j<N; j++){
			y[i] += Qt[i][j]*b[j];
		}
	}

	makeSolve(R, x, y, N);

	print_vector("X", x, N);

	for(i=0; i<N; i++){
		y[i] = 0;
		for(j=0; j<N; j++){
			y[i] += A[i][j]*x[j];
		}
	}

	print_vector("test", y, N);


	for(i=0; i<N; i++){
		free(A[i]);
		free(Q[i]);
		free(Qt[i]);
		free(R[i]);
	}
	free(A);
	free(Q);
	free(Qt);
	free(R);
	free(b);
	free(y);
	free(x);
	fclose(f);
}


void print_matrix(char *name, float **A, int m, int n){
	int i, j;

	printf("----%s----\n", name);
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}

void print_vector(char *name, float *A, int n){
	int j;

	printf("----%s----\n", name);
	for(j=0; j<n; j++){
		printf("%f ", A[j]);
	}
	printf("\n\n");

}

void makeHouseholder(float **A, float **Q, int N){
	int i, j, k;
	float norma;
	float *V;
	float **E;
	float **H;
	float **QQ;

	E  = (float **)malloc(N*sizeof(float *));
	H  = (float **)malloc(N*sizeof(float *));
	QQ = (float **)malloc(N*sizeof(float *));

	for(i=0 ; i<N; i++){
		E[i]  = (float *)malloc(N * sizeof(float));
		H[i]  = (float *)malloc(N * sizeof(float));
		QQ[i] = (float *)malloc(N * sizeof(float));
	}

	for(j=0; j<N; j++){
		for(k=0; k<N; k++){
			if(j==k){
				E[j][k] = 1.;
				Q[j][k] = 1.;
			}
			else{
				E[j][k] = 0.;
				Q[j][k] = 0.;
			}
		}
	}

	V = (float *)malloc(N*sizeof(float));

	for(i=0; i<N-1; i++){
//	makeV(i)
		for(j=0; j<i; j++) V[j] = 0.;

		V[i] = A[i][i] + sign(A[i][i])*makeNorma(A, i, N);

		for(j=i+1; j<N; j++) V[j] = A[j][i];

//  makeH
		norma = 0;
		for(j=0; j<N; j++) norma += V[j]*V[j];

		for(j=0; j<N; j++){
			for(k=0; k<N; k++){
				H[j][k] = E[j][k] - 2.*V[j]*V[k]/norma;
			}
		}

//  Q = Q*H
		mult_matrix(Q, H, QQ, N);
		for(j=0; j<N; j++){
			for(k=0; k<N; k++){
				Q[j][k] = QQ[j][k];
			}
		}

//  Ai+1 = H*Ai
		mult_matrix(H, A, QQ, N);
		for(j=0; j<N; j++){
			for(k=0; k<N; k++){
				A[j][k] = QQ[j][k];
			}
		}
	}

	for(i=0; i<N; i++){
		free(E[i]);
		free(H[i]);
		free(QQ[i]);
	}
	free(E);
	free(H);
	free(QQ);

	free(V);
}

float sign(float a){
	if (a==0)	return 0.;
	if (a>0)	return 1.;
	else 		return -1.;
}

float makeNorma(float **A, int k, int N){
	int i;
	float norma;

	norma = 0.0;
	for(i=k; i<N; i++){
		norma += A[i][k]*A[i][k];
	}

	return sqrt(norma);
}

void mult_matrix(float **A, float **B, float **C, int N){
	int i;
	int j;
	int k;

	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			C[i][j] = 0;
			for(k=0; k<N; k++){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void makeSolve(float **R, float *x, float *y, int N){
	int i, j;
	float sum;

	for(i=N-1; i>=0; i--){
		sum = 0;
		for(j=i+1; j<N; j++){
			sum += R[i][j]*x[j];
		}
		x[i] = (y[i] - sum)/R[i][i];
	}

}
