/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "problem_region_jacobian_residual_m4_p.h"
#include "main.h"
#include "problem_region_jacobian_residual_m4_p_terminate.h"
#include "problem_region_jacobian_residual_m4_p_initialize.h"

/* Function Declarations */
static void argInit_20x3x3_real_T(double result[180]);
static double argInit_real_T(void);
static void main_problem_region_jacobian_residual_m4_p(void);

/* Function Definitions */
static void argInit_20x3x3_real_T(double result[180])
{
  int idx0;
  int idx1;
  int idx2;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 20; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      for (idx2 = 0; idx2 < 3; idx2++) {
        /* Set the value of the array element.
           Change this value to the value that the application requires. */
        result[(idx0 + 20 * idx1) + 60 * idx2] = argInit_real_T();
      }
    }
  }
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_problem_region_jacobian_residual_m4_p(void)
{
  double dv0[180];
  double dv1[180];
  double solution[180];
  static double Jacobian[32400];
  double residual[180];

  /* Initialize function 'problem_region_jacobian_residual_m4_p' input arguments. */
  /* Initialize function input argument 'q'. */
  /* Initialize function input argument 'q_past'. */
  /* Call the entry-point 'problem_region_jacobian_residual_m4_p'. */
  argInit_20x3x3_real_T(dv0);
  argInit_20x3x3_real_T(dv1);
  problem_region_jacobian_residual_m4_p(dv0, dv1, argInit_real_T(),
    argInit_real_T(), solution, Jacobian, residual);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  problem_region_jacobian_residual_m4_p_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_problem_region_jacobian_residual_m4_p();

  /* Terminate the application.
     You do not need to do this more than one time. */
  problem_region_jacobian_residual_m4_p_terminate();
  return 0;
}

/* End of code generation (main.c) */
