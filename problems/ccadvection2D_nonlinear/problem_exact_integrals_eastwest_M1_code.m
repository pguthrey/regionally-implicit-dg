function [Jac_cell] = problem_exact_integrals_eastwest_M1_code(q,q_past,nuv1,nuf) %#codegen

Jac_cell = zeros(1,1);

u_1 = q(1);
u_past_1 = q_past(1);
Jac_cell(1,1) = u_1;
Jac_cell(1,1) = 4.0;
