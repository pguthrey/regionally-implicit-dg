function [ qplot ] = problem_cons2plot(qcons,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

e = qcons(4,:)/qcons(1,:);
v = qcons(3,:)/qcons(1,:);
u = qcons(2,:)/qcons(1,:);
s = sqrt(u^2 + v^2);
rho = qcons(1,:);
E = qcons(4,:);
p = (data.appdata.gamma-1)*(E-rho*s^2/2) ;
my = qcons(3,:);
mx = qcons(2,:);

qplot = [rho;mx;my;E;p;u;v;e;s];

end

