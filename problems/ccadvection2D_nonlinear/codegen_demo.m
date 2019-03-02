clc
clear all 
close all
mex -setup cpp
codegen problem_exact_integrals_eastwest_M1_code.m -args {zeros(1,1),zeros(1,1),1.0,1.0}
