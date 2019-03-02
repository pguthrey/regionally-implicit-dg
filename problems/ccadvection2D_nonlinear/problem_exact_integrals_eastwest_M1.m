function [trunceast,truncwest,residual_cell, ...
           Jac_trunc_east,Jac_trunc_west,Jac_cell] = problem_exact_integrals_eastwest_M1(q,q_past,nuv1,nuf) 
u_1 = q(1);
u_past_1 = q_past(1);
trunceast(1,1) = 0;
trunceast(1) = nuv1*(4.0*nuf*u_1);
truncwest(1,1) = 0;
truncwest(1) = nuv1*(-4.0*nuf*u_1);
residual_cell(1,1) = 0;
residual_cell(1) = 4.0*u_1 - 4.0*u_past_1;
Jac_trunc_east(1,1) = 0;
Jac_trunc_east(1,1) = nuv1*(4.0*nuf);
Jac_trunc_west(1,1) = 0;
Jac_trunc_west(1,1) = nuv1*(-4.0*nuf);
Jac_cell(1,1) = 0;
Jac_cell(1,1) = 4.0;
