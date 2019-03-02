/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * problem_region_Jacobian_residual_M4_P.h
 *
 * Code generation for function 'problem_region_Jacobian_residual_M4_P'
 *
 */

#ifndef PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_H
#define PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "problem_region_Jacobian_residual_M4_P_types.h"

/* Function Declarations */
extern void problem_region_Jacobian_residual_M4_P(const double q[30], const
  double q_past[30], double nuv1, double Jacobian[900], double residual[30]);

#endif

/* End of code generation (problem_region_Jacobian_residual_M4_P.h) */
