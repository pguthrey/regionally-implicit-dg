function [DGsolution,auxiliary,data] = problem_initialize_auxiliary(DGsolution,data)

auxiliary.pressure_hartree = zeros(data.Ls,data.Nv1,data.Nv2,data.Nv3);

[data] = initialize_gauss_hermite_quadrature(data,data.M);