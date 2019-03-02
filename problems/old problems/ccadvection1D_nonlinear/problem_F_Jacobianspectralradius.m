function [ rho ] = problem_F_Jacobianspectralradius(q,quadpoint,appdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Linear Advection
rho = abs(appdata.nuf);

end

