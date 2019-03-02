function [ rho ] = problem_G_Jacobianspectralradius(q,quadpoint,appdata)

v1 = quadpoint(1);
v2 = quadpoint(2);
Rm = 1;
Rq = 1;
gamma = 1;%sqrt(Rm^2-v1.^2-v2.^2);
rho = 0;%Rq*abs(appdata.E2-v1.*appdata.B3./gamma);
end

