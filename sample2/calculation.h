#ifndef __CALCULATION_H_INCLUDED__
#define __CALCULATION_H_INCLUDED__

typedef struct{
	int size;
	double *v;
}Vector;

typedef struct{
	int height;
	int width;
	double *a;
}Matrix;

Vector *CreateVector(int size);
void FreeVector(Vector* vector);
double GetDotProduct(Vector* vector1, Vector* vector2);

Matrix *CreateMatrix(int width, int height);
void FreeMatrix(Matrix* matrix);

//$B%Y%/%H%k(Bx$B$rBh0lMWAG$,(Bx$B$N%N%k%`(B(sign(x[0]))$B$G;D$j$,(B0$B$N(Bx'$B%Y%/%H%k$KJQ49(B
//x'[0]$B$rJV$9!#(B
double Householder(Vector *x);
//$B9TNs(Bm$B$r;0=EBP3Q2=$9$k!#(B
void Tridiag(Matrix *m);
//QR$BK!$K$h$C$FBP3Q2=!#BP3Q@.J,$K8GM-CM!#8m:9$O(Beps$B0JFb(B
void QR(Matrix* m, double eps);
#endif /* __CALCULATION_H_INCLUDED__ */
