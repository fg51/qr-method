#include "householder.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define EPS (1.0E-8)


bool check_calloc_error(double u[], double d[], double ds[]);

static void init_u(double *, const double *, const int k, const int n);
static double cal_sum(const double xs[], const int begin, const int num);
static bool create_u(double u[], const double xs[], const int k, const int num);
static void init_d(double d[], const double xs[], const double u[], const int k, const int length);
static void init_ds(double ds[], const double xs[], const double u[], const int k, const int length);
static double product_sum(const double xs[], const double ys[], const int begin, const int end);
static void update_d(double ys[], const double xs[], const int begin, const int end);
static void update_hessenberg(double [], const double [], const double [], const double [], const int);


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

