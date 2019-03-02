
DGregion = zeros(56,3,3);
DGregion_past = zeros(56,3,3);
nuv1 = .7;
nuv2 = .7;
[Jacobian,residual] = problem_region_Jacobian_residual_M6_P(DGregion,DGregion_past,nuv1,nuv2); 
