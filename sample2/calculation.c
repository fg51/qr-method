#include<stdio.h>
#include<math.h>
#include"calculation.h"

//Tridiagで使う。householderによってできたベクトルによる直交行列を右から掛ける
static void Mult(Matrix* x, Vector* v);
//行列mを単位行列に
static void InitMatrix(Matrix* m);
//2×2行列の固有値を取得し、a22成分に近い方を返す
static double GetEigenvalue22(double a11, double a12, double a21, double a22);
//QR法で使う。Aの左からP、Qの右からPを掛ける
static void NextAQ(Matrix* a, Matrix* q, int nowsize);

Vector *CreateVector(int size)
{
	Vector *temp;

	if((temp = (Vector*)malloc(sizeof(Vector))) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	if((temp->v = (double*)malloc(sizeof(double)*size)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	temp->size = size;

	return temp;
}

void FreeVector(Vector* vector)
{
	free(vector->v);
	free(vector);
}

double GetDotProduct(Vector* vector1, Vector* vector2)
{
	double sum;
	int i;

	if(vector1->size != vector2->size){
		fprintf(stderr, "Cannot get dot product.\n");
		exit(1);
	}

	for(i=0, sum=.0; i<vector1->size; i++){
		sum += vector1->v[i] * vector2->v[i];
	}

	return sum;
}

Matrix *CreateMatrix(int width, int height)
{
	Matrix *temp;

	if((temp = (Matrix*)malloc(sizeof(Matrix))) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	if((temp->a = (double*)malloc(sizeof(double)*width*height)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	temp->width = width;
	temp->height = height;

	return temp;
}

void FreeMatrix(Matrix* matrix)
{
	free(matrix->a);
	free(matrix);
}

double Householder(Vector *x)
{
	double norm;
	double weight;
	int i, j;
	

	norm = sqrt(GetDotProduct(x, x));

	if(norm){
		if(x->v[0] < 0) norm = - norm;
		x->v[0] += norm;

		weight = 1/(sqrt(2*norm*x->v[0]));

		for(i=0; i<x->size; i++){
			x->v[i] *= weight;
		}
	}

	return -norm;
}


void Tridiag(Matrix *m)
{
	Vector *v;
	int i, j;
	double tmp;

	v = CreateVector(m->width);

	for(i=0; i<m->width-2; i++){
		v->size--;
		for(j=i+1; j<m->width; j++){
			v->v[j-i-1] = m->a[i*m->width + j];
		}
	  if(!(tmp = Householder(v))) continue;

		Mult(m, v);

		m->a[i*m->width + i+1] = m->a[(i+1)*m->width + i] = tmp;
		for(j=i+2; j<m->width; j++){
			m->a[i*m->width + j] = m->a[j*m->width + i] = 0;
		}
		

	}

	FreeVector(v);
}

static void Mult(Matrix* x, Vector* v)
{
	int i, j, offset;
	double tmp;
	Vector *g;
  
	offset = x->width - v->size;
	g = CreateVector(v->size);

	for(i=offset; i<x->width; i++){
		g->v[i-offset] = 0;
		for(j=offset; j<x->width; j++){
			g->v[i-offset] += x->a[i*x->width + j]*v->v[j-offset];
		}
	}

	tmp = GetDotProduct(g, v);

	for(j=0; j<g->size; j++){
		g->v[j] = 2*(g->v[j] - v->v[j]*tmp);
	}
	for(i=offset; i<x->width; i++){
		for(j=offset; j<x->width; j++){
			x->a[i*x->width +j] -= (v->v[i-offset]*g->v[j-offset] + g->v[i-offset]*v->v[j-offset]);
		}
	}

	FreeVector(g);
}

void QR(Matrix* m, double eps)
{
	int nowsize;
	double eig1, eig2;
	int i, j;
	Matrix* q;
	Vector* v;

	Tridiag(m);

	nowsize = m->width;

	q = CreateMatrix(nowsize, nowsize);
	v = CreateVector(nowsize);

	while(nowsize > 1){
		double det, b, u;

		if(fabs(m->a[(nowsize-1)*m->width + nowsize-2]) < eps){
			nowsize--;
			continue;
		}

		u = GetEigenvalue22(m->a[(nowsize-2)*m->width + nowsize-2], m->a[(nowsize-2)*m->width + nowsize-1], 
												m->a[(nowsize-1)*m->width + nowsize-2], m->a[(nowsize-1)*m->width + nowsize-1]);

		for(i=0; i<nowsize; i++){
			m->a[i*m->width + i] -= u;
		}

		q->height = nowsize;
		q->width = nowsize;

		InitMatrix(q);

		NextAQ(m, q, nowsize);

		for(i=0; i<nowsize; i++){
			double sum=0;
			int k;
			
			for(j=0; j<nowsize; j++){
				for(k=i; k<nowsize; k++){
					sum += m->a[i*m->width + k]*q->a[k*q->width + j];
				}
				v->v[j] = sum;
				sum = 0;
			}
			for(j=0; j<nowsize; j++){
				m->a[i*m->width + j] = v->v[j];
			}
		}

		for(i=0; i<nowsize; i++){
			m->a[i*m->width + i] += u;
		}
	
	}

	FreeMatrix(q);
	FreeVector(v);
}

static void InitMatrix(Matrix* m)
{
	int i, j;

	for(i=0; i<m->height; i++){
		for(j=0; j<m->width; j++){
			if(i==j) m->a[i*m->width + j] = 1;
			else m->a[i*m->width + j] = 0;
		}
	}
}

static double GetEigenvalue22(double a11, double a12, double a21, double a22)
{
	double b, det, eig1, eig2, u;

	b = a22 + a11;
	det = b*b - 4*(a22*a11 - a12*a21);

	if(det < 0) det = 0;
	eig1 = (b + sqrt(det))/2;
	eig2 = (b - sqrt(det))/2;

	if(fabs(a22 - eig1) < fabs(a22 - eig2)){
		u = eig1;
	}else{
		u = eig2;
	}

	return u;
}

static void NextAQ(Matrix* a, Matrix* q, int nowsize)
{
	int i, j;
	double alpha, s, c;
	double temp;

	for(i=0; i<nowsize-1; i++){
		alpha = sqrt(a->a[i*a->width + i]*a->a[i*a->width + i] + a->a[(i+1)*a->width + i]*a->a[(i+1)*a->width + i]);
		s = alpha ? a->a[(i+1)*a->width + i]/alpha : alpha;
		c = alpha ? a->a[i*a->width + i]/alpha : alpha;
		
		for(j=i+1; j<nowsize;  j++){
			temp = -a->a[i*a->width + j]*s + a->a[(i+1)*a->width + j]*c;
			a->a[i*a->width + j] = a->a[i*a->width + j]*c + a->a[(i+1)*a->width + j]*s;
			a->a[(i+1)*a->width + j] = temp;
		}

		for(j=0; j<nowsize; j++){
			temp = -q->a[j*q->width + i]*s + q->a[j*q->width + i+1]*c;
			q->a[j*q->width + i] = q->a[j*q->width + i]*c + q->a[j*q->width + i+1]*s;
			q->a[j*q->width + i+1] = temp;
		}
		a->a[i*a->width + i] = alpha;
		a->a[(i+1)*a->width + i] = 0;
	}
}

