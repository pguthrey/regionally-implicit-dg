/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.cpp
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
#include "problem_exact_integrals_eastwest_M1_code.h"
#include "main.h"
#include "problem_exact_integrals_eastwest_M1_code_terminate.h"
#include "problem_exact_integrals_eastwest_M1_code_initialize.h"

/* Function Declarations */
static double argInit_real_T();
static void main_problem_exact_integrals_eastwest_M1_code();

/* Function Definitions */
static double argInit_real_T()
{
  return 0.0;
}

static void main_problem_exact_integrals_eastwest_M1_code()
{
  double Jac_cell;

  /* Initialize function 'problem_exact_integrals_eastwest_M1_code' input arguments. */
  /* Call the entry-point 'problem_exact_integrals_eastwest_M1_code'. */
  Jac_cell = problem_exact_integrals_eastwest_M1_code(argInit_real_T(),
    argInit_real_T(), argInit_real_T(), argInit_real_T());
}

int main(int, const char * const [])
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  problem_exact_integrals_eastwest_M1_code_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_problem_exact_integrals_eastwest_M1_code();

  /* Terminate the application.
     You do not need to do this more than one time. */
  problem_exact_integrals_eastwest_M1_code_terminate();
  return 0;
}

/* End of code generation (main.cpp) */
