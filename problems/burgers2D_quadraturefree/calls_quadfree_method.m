
DGregion = zeros(20,3,3);
DGregion_past = zeros(20,3,3);
nuv1 = .7;
nuv2 = .7;
[solution,Jacobian,residual] = problem_region_jacobian_residual_m4_p(DGregion,DGregion_past,nuv1,nuv2); 
