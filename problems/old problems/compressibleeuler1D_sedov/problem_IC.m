function [ qIC ] = problem_IC(point,appdata)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------re

x = point(1);
%NonCC Advection
%%{
Pref =  .1;
Pcore = 10;
radius = .3;
delta = 1/100;
p = Pref + (Pcore-Pref)/2*(1+tanh((radius/2+x)/delta).*tanh((radius/2-x)/delta));

rho = 1+0.*x;
u = 0.*x;
E = p/(appdata.gamma-1)+rho*u^2/2;
qIC(3,1) = E;
qIC(2,1) = rho*u;
qIC(1,1) = rho;



end

