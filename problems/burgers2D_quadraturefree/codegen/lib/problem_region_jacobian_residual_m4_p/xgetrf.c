/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgetrf.c
 *
 * Code generation for function 'xgetrf'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "problem_region_jacobian_residual_m4_p.h"
#include "xgetrf.h"
#include "xswap.h"

/* Function Definitions */
void xgetrf(double A[32400], int ipiv[180], int *info)
{
  int i0;
  int j;
  int c;
  int jA;
  int ix;
  double smax;
  int jy;
  double s;
  int b_j;
  int ijA;
  for (i0 = 0; i0 < 180; i0++) {
    ipiv[i0] = 1 + i0;
  }

  *info = 0;
  for (j = 0; j < 179; j++) {
    c = j * 181;
    jA = 1;
    ix = c;
    smax = fabs(A[c]);
    for (jy = 2; jy <= 180 - j; jy++) {
      ix++;
      s = fabs(A[ix]);
      if (s > smax) {
        jA = jy;
        smax = s;
      }
    }

    if (A[(c + jA) - 1] != 0.0) {
      if (jA - 1 != 0) {
        ipiv[j] = j + jA;
        xswap(A, j + 1, j + jA);
      }

      i0 = (c - j) + 180;
      for (jA = c + 1; jA < i0; jA++) {
        A[jA] /= A[c];
      }
    } else {
      *info = j + 1;
    }

    jA = c;
    jy = c + 180;
    for (b_j = 1; b_j <= 179 - j; b_j++) {
      smax = A[jy];
      if (A[jy] != 0.0) {
        ix = c + 1;
        i0 = (jA - j) + 360;
        for (ijA = 181 + jA; ijA < i0; ijA++) {
          A[ijA] += A[ix] * -smax;
          ix++;
        }
      }

      jy += 180;
      jA += 180;
    }
  }

  if ((*info == 0) && (!(A[32399] != 0.0))) {
    *info = 180;
  }
}

/* End of code generation (xgetrf.c) */
