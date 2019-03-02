function [Jac_cell,residual_cell] = problem_exact_jaccell_M1_P(q,q_past,nuv1,nuv2) 
u_1 = q(1);
u_past_1 = q_past(1);
residual_cell(1,1) = 0;
residual_cell(1,1) = 4.0*u_1 - 4.0*u_past_1;
Jac_cell(1,1) = 0;
Jac_cell(1,1) = 4.0;
