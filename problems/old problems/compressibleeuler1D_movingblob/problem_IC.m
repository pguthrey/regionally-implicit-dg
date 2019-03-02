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
rho = exp(-200.0*(x-.5)^2)+.1;
u = appdata.fluidu;
gamma = appdata.gamma;
p = appdata.fluidp;

E = p/(gamma-1) + rho*u^2/2;
qIC(3,1) = E;
qIC(2,1) = rho*u;
qIC(1,1) = rho;

end

