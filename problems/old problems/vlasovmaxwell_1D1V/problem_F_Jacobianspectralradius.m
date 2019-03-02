function [ rho ] = problem_F_Jacobianspectralradius(q,quadpoint,appdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Rq = -1;
v1 = quadpoint(1);
gamma = 1;%sqrt(Rm^2-v1.^2-v2.^2);
check = (v1 < appdata.v1_ub)*(v1 > appdata.v1_lb);
rho = abs(Rq*appdata.E1)*check;

end

