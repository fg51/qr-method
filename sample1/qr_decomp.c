#include "qr_decomp.h"

#include <stdlib.h>
#include <math.h>


#define EPS (1.0E-8)


static double cal_wa(const double x00, const double x01, const double x10, const double x11);
static double cal_mu(const double a00, const double a01, const double a10, const double a11);
static void init_q(double q[], const int m);
static void update_a(double a[], const double sinx, const double cosx, const int begin, const int end, const int n);
static void update_q(double q[], const int end, const int i, const double sinx, const double cosx);
static double cal_sum_w_q(const double w[], const double q[], const int begin, const int end, const int j);


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
