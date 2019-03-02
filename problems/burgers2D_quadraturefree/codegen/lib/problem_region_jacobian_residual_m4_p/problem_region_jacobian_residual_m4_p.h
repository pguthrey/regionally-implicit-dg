/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * problem_region_jacobian_residual_m4_p.h
 *
 * Code generation for function 'problem_region_jacobian_residual_m4_p'
 *
 */

#ifndef PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_H
#define PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "problem_region_jacobian_residual_m4_p_types.h"

/* Function Declarations */
extern void problem_region_jacobian_residual_m4_p(const double q[180], const
  double q_past[180], double nuv1, double nuv2, double solution[180], double
  Jacobian[32400], double residual[180]);

#endif

/* End of code generation (problem_region_jacobian_residual_m4_p.h) */
