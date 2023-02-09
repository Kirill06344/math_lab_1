#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cmath.h"

double function(double x)
{
  return 1 / (1 + x);
}


void fill_values(double* x_values, double* y_values, int nodes)
{
  for (int k = 0; k < nodes; ++k) {
    x_values[k] = 0.1 * k;
    y_values[k] = function(x_values[k]);
  }
}

double interpolateLagrange(double x, double* x_values, double* y_values, int nodes)
{
  double res = 0;
  for (int i = 0; i < nodes; ++i) {
    double term = 1;
    for (int k = 0; k < nodes; ++k) {
      if (k == i) continue;
      term *= (x - x_values[k]) / (x_values[i] - x_values[k]);
    }
    res += term * y_values[i];
  }
  return res;
}

void printResults(double* x_values, double* y_values, double* b, double* c, double* d, int nodes)
{
  int last;
  printf("Results...\n");
  printf(" x       f(x)    Lagrange   spline\n");
  printf("----------------------------------\n");
  for (int k = 0; k < 10; ++k) {
    double x = 0.05 + 0.1 * k;
    double spline = seval(nodes, x, x_values, y_values, b, c, d, &last);
    printf("%lf %lf %lf %lf\n", x, function(x), interpolateLagrange(x, x_values, y_values, nodes), spline);
  }
}


double f1 (double x)
{
  return pow(fabs(x - tan(x)), -1);
}

double f2 (double x)
{
  return pow(fabs(x - tan(x)), -0.5);
}


void quanc8_results(double (*f)(double x), double l, double r, double epsabs, double epsrel) {
  double result, errest, posn, flag;
  int nfe;

  quanc8(f,l, r, epsabs, epsrel, &result, &errest, &nfe, &posn, &flag);
  printf ("%s\n\n", cmathmsg(QUANC8_C, flag));

  printf ("nfe = %d \n", nfe);
  printf ("integral = %lf\n", result);
  printf ("error = %lf \n", errest);
  printf("flag = %lf\n", flag);
  if (flag < 0)
  {
    printf ("trouble spot at x = %e \n", posn);
    printf ("%d unconverged subintervals\n", abs(flag));
  }
}

int main()
{
  #define ndim 20
  const int nodes = 11;
  double x_values[ndim], y_values[ndim];
  double b[ndim], c[ndim], d[ndim];
  int flag;

  fill_values(x_values, y_values, nodes);


  //Using natural end condition at x without slope.That's why arguments e1,e2,s1,s2 = 0.
  spline(nodes, 0, 0, 0, 0, x_values, y_values, b, c, d, &flag);
  printf("%s\n\n", cmathmsg(SPLINE_C, flag));

  if (flag == 0) {
    printResults(x_values, y_values, b, c, d, nodes);
  }


  double l = 2.00;
  double r = 5.00;
  double epsrel = 0.0;
  double epsabs = 1.0e-3;

  quanc8_results(f1, l, r, epsabs, epsrel);
  quanc8_results(f2, l, r, epsabs, epsrel);



}

