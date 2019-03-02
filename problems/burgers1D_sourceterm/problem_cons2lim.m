function [lim_pos,lim_max,lim_min] = problem_cons2lim(qcons,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

error('whoops')

e = qcons(3,:)/qcons(1,:);
u = qcons(2,:)/qcons(1,:);
rho = qcons(1,:);
E = qcons(3,:);
p = (data.appdata.gamma-1)*(E-rho*u^2/2) ;
mx = qcons(2,:);

lim_pos = NaN*[rho;mx;E];
lim_max = [rho;mx;E];
lim_min = [rho;mx;E];

%{
lim_pos = [rho;E;p;e];
lim_max = [E;p;e];%[rho;u;v;E;p;mx;my;];
lim_min = [E;p;e];%lim_max;
%}
end

