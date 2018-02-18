/*
 * QR method
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


#define EPS (1.0E-8)
#define N (4)


bool check_calloc_error(double u[], double d[], double ds[]);


bool householder(double *, const int);
static void init_u(double *, const double *, const int k, const int n);
static double cal_sum(const double xs[], const int begin, const int num);
static bool create_u(double u[], const double xs[], const int k, const int num);
static void init_d(double d[], const double xs[], const double u[], const int k, const int length);
static void init_ds(double ds[], const double xs[], const double u[], const int k, const int length);
static double product_sum(const double xs[], const double ys[], const int begin, const int end);
static void update_d(double ys[], const double xs[], const int begin, const int end);
static void update_hessenberg(double [], const double [], const double [], const double [], const int);

void qr_method(double *, int);
static double cal_wa(const double x00, const double x01, const double x10, const double x11);
static double cal_mu(const double a00, const double a01, const double a10, const double a11);
static void init_q(double q[], const int m);
static void update_a(double a[], const double sinx, const double cosx, const int begin, const int end, const int n);
static void update_q(double q[], const int end, const int i, const double sinx, const double cosx);
static double cal_sum_w_q(const double w[], const double q[], const int begin, const int end, const int j);

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


/*
 * to upper hessenberg matrix with householder
 */
bool householder(double xs[], const int num){
  double *u = (double *)calloc(num, sizeof(double));
  double *d = (double *)calloc(num, sizeof(double));
  double *ds =(double *)calloc(num, sizeof(double));
  if (check_calloc_error(u, d, ds) == true) {
    return true;
  }

  for (int k = 0; k <= num - 3; k++) {
    //  変換行列 H の構築
    if (create_u(u, xs, k, num) == true) {
      continue;
    }

    //  similarity transform
    init_d(d, xs, u, k, num);
    update_d(d, u, k + 1, num);

    init_ds(ds, xs, u, k, num);
    update_d(ds, u, k + 1, num);

    update_hessenberg(xs, u, d, ds, num);
  }
  return false;
}

bool check_calloc_error(double u[], double d[], double ds[])
{
  if (u == NULL) {
    puts("u is calloc error");
    return true;
  }
  if (d == NULL) {
    puts("d is calloc error");
    return true;
  }
  if (ds == NULL) {
    puts("ds is calloc error");
    return true;
  }
  return false;
}


static void init_u(double ys[], const double xs[], const int k, const int num)
{
 for (int i = 0; i <= k; i++) {
   ys[i] = 0.0;
 }
 for (int i = k + 1; i < num; i++) {
   ys[i] = xs[num * i + k];
 }
}


static bool create_u(double u[], const double xs[], const int k, const int num)
{
  init_u(u, xs, k, num);

  double sum = cal_sum(u, k, num);
  if (fabs(u[k+1]) < EPS) {
    return true;
  }

  const double sigma = sqrt(sum) * u[k + 1] / fabs(u[k + 1]);
  u[k+1] += sigma;
  const double v_norm = sqrt(2.0 * sigma * u[k + 1]);
  for (int i = k + 1; i < num; i++) {
    u[i] /= v_norm;
  }
  return false;
}


static double cal_sum(const double xs[], const int begin, const int num)
{
  double sum = 0.0;
  for (int i = begin + 1; i < num; i++) {
    sum += xs[i] * xs[i];
  }
  return sum;
}


static void init_d(double d[], const double xs[], const double u[], const int k, const int length)
{
  for (int i = 0; i < length; i++) {
    d[i] = 0.0;
    for (int j = k + 1; j <= length - 1; j++) {
      d[i] += xs[length * i + j] * u[j];
    }
  }
}


static void init_ds(double ds[], const double xs[], const double u[], const int k, const int length)
{
  for (int i = 0; i < length; i++) {
    ds[i] = 0.0;
    for (int j = k + 1; j <= length - 1; j++) {
      ds[i] += xs[length * j + i] * u[j];
    }
  }
}


static double product_sum(const double xs[], const double ys[], const int begin, const int end)
{
  double total = 0.0;
  // for (int i = k + 1; i < end; i++) {
  for (int i = begin; i < end; i++) {
     total += xs[i] * ys[i];
  }
  return total;
}

static void update_d(double ys[], const double xs[], const int begin, const int end)
{
  const double ud  = product_sum(xs, ys,  begin, end);
  for (int i = 0; i < end; i++) {
     ys[i] = 2.0 * (ys[i] - ud * xs[i]);
  }
}

static void update_hessenberg(double xs[], const double u[], const double d[], const double ds[], const int num)
{
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < num; j++) {
      xs[num * i + j] -= u[i] * ds[j] + d[i] * u[j];
    }
  }
}


void qr_method(double a[],int n){
  double *q = (double *)calloc(n * n, sizeof(double));
  double *w = (double *)calloc(n,     sizeof(double));

  int m = n;

  while ( m != 1 ) {
    if (fabs(a[n * (m - 1) + m - 2]) < EPS) {
      m -= 1;
      continue;
    }

    // 原点移動 mu
    const double mu = cal_mu(
      a[n * (m - 2) + m - 2], a[n * (m - 2) + m - 1],
      a[n * (m - 1) + m - 2], a[n * (m - 1) + m - 1]);
    for (int i = 0; i < m; i++) {
      a[n * i + i] -= mu;
    }

    // QR method
    init_q(q, m);

    for (int i = 0; i < m-1; i++) {
      const double sum1 = sqrt(a[n * i + i] * a[n * i + i] + a[n * i + n + i] * a[n * i + n + i]);
      const double sinx = (fabs(sum1) < EPS) ? 0.0 : a[n * i + n + i] / sum1;
      const double cosx = (fabs(sum1) < EPS) ? 0.0 : a[n * i + i] / sum1;

      update_a(a, sinx, cosx, i, m, n);
      a[n * i + i] = sum1;
      update_q(q, m, i, sinx, cosx);
    }

    for (int i = 0; i < m; i++) {
      for(int j = i; j < m; j ++) {
        w[j] = a[n * i + j];
      }
      for (int j = 0; j < m; j++) {
        a[n * i + j] = cal_sum_w_q(w, q, i, m, j);
      }
    }
    for(int i = 0; i < m; i++) {
      a[n * i + i] += mu;
    }
  }
}


static double cal_wa(const double x00, const double x01, const double x10, const double x11)
{
  const double sum1 = x00 + x11;
  const double sum2 = x00 * x11 - x01 * x10;
  const double wa = sum1 * sum1 - 4.0 * sum2;
  return (wa < 0.0) ? 0.0 : sqrt(wa);
}


static double cal_mu(const double a00, const double a01, const double a10, const double a11)
{
    const double wa = cal_wa(a00, a01, a10, a11);

    const double sum1 = a00 + a11;
    const double sum2 = a00 * a11 - a01 * a10;
    const double lam1 = 0.5 * (sum1 + wa);
    const double lam2 = sum2 / lam1;
    return (fabs(a11 - lam1) < fabs(a11 - lam2)) ? a11 - lam1 : a11 - lam2;
}


static void init_q(double q[], const int m)
{
  for(int i = 0; i < m * m; i++) {
    q[i] = 0.0;
  }
  for(int i = 0; i < m; i++) {
    q[m * i + i] = 1.0;
  }
}


static void update_a(double a[], const double sinx, const double cosx, const int begin, const int end, const int n)
{
  for (int j = begin + 1; j < end; j++) {
    const int nij = n * begin + j;
    const double sum2 = a[nij] * cosx + a[nij + n] * sinx;
    a[nij + n] = -a[nij] * sinx + a[nij + n] * cosx;
    a[nij] = sum2;
  }
  a[n * begin + n + begin] = 0.0;
}


static void update_q(double q[], const int end, const int i, const double sinx, const double cosx)
{
  for (int j = 0; j < end; j++) {
    const double sum2 = q[end * j + i] * cosx + q[end * j + i + 1] * sinx;
    q[end * j + i + 1] = -q[end * j + i] * sinx + q[end * j + i + 1] * cosx;
    q[end * j + i] = sum2;
  }
}


static double cal_sum_w_q(const double w[], const double q[], const int begin, const int end, const int j)
{
  double total = 0.0;
  for (int k = begin; k < end; k++) {
    total += w[k] * q[end * k + j];
  }
  return total;
}
