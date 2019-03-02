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
#include "problem_region_Jacobian_residual_M4_P.h"
#include "main.h"
#include "problem_region_Jacobian_residual_M4_P_terminate.h"
#include "problem_region_Jacobian_residual_M4_P_initialize.h"

/* Function Declarations */
static void argInit_10x3_real_T(double result[30]);
static double argInit_real_T(void);
static void main_problem_region_Jacobian_residual_M4_P(void);

/* Function Definitions */
static void argInit_10x3_real_T(double result[30])
{
  int idx0;
  double result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 10; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_real_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 10] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 20] = argInit_real_T();
  }
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_problem_region_Jacobian_residual_M4_P(void)
{
  double q_tmp[30];
  double Jacobian[900];
  double residual[30];

  /* Initialize function 'problem_region_Jacobian_residual_M4_P' input arguments. */
  /* Initialize function input argument 'q'. */
  argInit_10x3_real_T(q_tmp);

  /* Initialize function input argument 'q_past'. */
  /* Call the entry-point 'problem_region_Jacobian_residual_M4_P'. */
  problem_region_Jacobian_residual_M4_P(q_tmp, q_tmp, argInit_real_T(), Jacobian,
    residual);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  problem_region_Jacobian_residual_M4_P_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_problem_region_Jacobian_residual_M4_P();

  /* Terminate the application.
     You do not need to do this more than one time. */
  problem_region_Jacobian_residual_M4_P_terminate();
  return 0;
}

/* End of code generation (main.c) */
