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

//ベクトルxを第一要素がxのノルム(sign(x[0]))で残りが0のx'ベクトルに変換
//x'[0]を返す。
double Householder(Vector *x);
//行列mを三重対角化する。
void Tridiag(Matrix *m);
//QR法によって対角化。対角成分に固有値。誤差はeps以内
void QR(Matrix* m, double eps);
#endif /* __CALCULATION_H_INCLUDED__ */
