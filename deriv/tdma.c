/*
         FILE: tdma.c
       AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
  DESCRIPTION: C version of TDMA
 */

void tdma(double *a, double *b, double *c, double *rhs, int n)
{
  int i;
  double m;

  // Forward elimination
  for(i = 1; i < n; i++)
  {
    m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    rhs[i] -= m * rhs[i - 1];
  }

  // Backward substitution
  rhs[n - 1] /= b[n - 1];
  for(i = n - 2; i >= 0; i--)
  {
    rhs[i] -= c[i] * rhs[i + 1];
    rhs[i] /= b[i];
  }
}
