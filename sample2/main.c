#include<stdio.h>
#include<math.h>
#include"calculation.h"

#define SIZE 3
#define EPS 0.0000001

void DisplayMatrix(Matrix *m);

int main()
{
	Matrix *m;
	int  j;

	m = CreateMatrix(SIZE, SIZE);

	m->a[0*m->width + 0] = 0.0;
	m->a[0*m->width + 1] = 2.0;
	m->a[0*m->width + 2] = 2.0;

	m->a[1*m->width + 0] = 2.0;
	m->a[1*m->width + 1] = 1.0;
	m->a[1*m->width + 2] = 0.0;

	m->a[2*m->width + 0] = 2.0;
	m->a[2*m->width + 1] = 0.0;
	m->a[2*m->width + 2] = -1.0;

	printf("Original Matrix:\n");
	DisplayMatrix(m);

	QR(m, EPS);

	printf("QR Diag Matrix:\n");
	DisplayMatrix(m);

	printf("QR EigenValue;\n");
	for(j=0; j<m->height; j++){
		printf("%f ", m->a[j*m->width + j]);
	}
	printf("\n\n");

	FreeMatrix(m);

	return 0;
}


void DisplayMatrix(Matrix *m)
{
	int i, j;

	for(i=0; i<m->height; i++){
		for(j=0; j<m->width; j++){
			printf("%f ", m->a[i*m->width + j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
}

