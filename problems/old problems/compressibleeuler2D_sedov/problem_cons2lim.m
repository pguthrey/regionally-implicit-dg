function [lim_pos,lim_max,lim_min,all] = problem_cons2lim(qcons,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


v = qcons(3,:)/qcons(1,:);
u = qcons(2,:)/qcons(1,:);
N = u^2 + v^2;
rho = qcons(1,:);
E = qcons(4,:);
p = (data.appdata.gamma-1)*(E- rho*N/2) ;
my = qcons(3,:);
mx = qcons(2,:);
e = p/(data.appdata.gamma-1)/rho;

all = [rho;u;v;E;p;mx;my;e];

u = NaN;
v = NaN;
mx = NaN;
my = NaN;

lim_pos = [rho;u;v;E;p;mx;my;e];

rho = NaN;
lim_max = [rho;u;v;E;p;mx;my;e];
lim_min = [rho;u;v;E;p;mx;my;e];


end

