function [ rho ] = problem_G_Jacobianspectralradius(q,quadpoint,appdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%CC Advection
%rho = 1;
%Non-CC Advection
gamma = sqrt(1+quadpoint(1)^2+quadpoint(2)^2);

rho = abs(quadpoint(1)/gamma);%abs(quadpoint(1));%*(max(abs(quadpoint(1:2))) < 1);
end

