function [ rho ] = problem_F_Jacobianspectralradius(q,quadpoint,appdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Linear Advection
gamma = sqrt(1+quadpoint(1)^2+quadpoint(2)^2);

rho = abs(quadpoint(2)/gamma);%abs(quadpoint(2));%*(max(abs(quadpoint(1:2))) < 1);

end

