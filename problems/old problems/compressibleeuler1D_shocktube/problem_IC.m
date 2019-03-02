function [ qIC ] = problem_IC(point,appdata)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------re

x = point(1);
y = point(2);
%NonCC Advection

%http://www.clawpack.org/_static/amrclaw/examples/euler_2d_quadrants/qinit.f.html

%%{

p = (x<.5)+(x>=.5)*.1;
rho = (x<.5)+(x>=.5)*.125;
u = 0.*x;
E = p/(appdata.gamma-1)+rho*u^2/2;
qIC(3,1) = E;
qIC(2,1) = rho*u;
qIC(1,1) = rho;

end

