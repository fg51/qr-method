/*
 * QR method
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


#define EPS (1.0E-8)
#define N (4)

#include "householder.h"
#include "qr_decomp.h"


void print_matrix(const double [], const int, const int);


int main(void) {
  // ===== 行列(入力データ) ================================================
  //   double a[N*N]={4.0,-6.0,5.0, -6.0,3.0,4.0, 5.0,4.0,-3.0}; 
  //   double a[N*N]={-1.,2.,3.,3., 2.,-3.,4.,1., 3.,4.,-1.,2., 3.,1.,2.,-3.};
  double xs[N * N] = {
    0.2, .32, .12, .3,
    0.1, .15, .24, .32,
    0.2, .24, .46, .36,
    0.6, .4,  .32, .2,
  };
  // =======================================================================
  printf(" matrix \n");
  print_matrix(xs, N, N);

  if (householder(xs, N)) {
    return 1;
  }
  printf(" after house holder \n");
  print_matrix(xs, N, N);

  qr_method(xs, N );
  printf(" after QR method（対角項が固有値） \n");
  print_matrix(xs, N, N);

  return 0;
}


/// 入力；   a= 行列, n= 表示する列数, m= 表示する行数
// void show_matrix( double a[], int column, int row)
void print_matrix(const double xs[], const int n, const int m)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("  %10.6f", xs[i * m + j]);
    }
    printf("\n");
  }
}



