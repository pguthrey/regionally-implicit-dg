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
p = .029*(x<.8)*(y<.8)  ...
    +.3*(x>=.8)*(y<.8)  ...
    +.3*(x<.8)*(y>=.8)  ...
    +1.5*(x>=.8)*(y>=.8);

rho = .138*(x<.8)*(y<.8) ...
    +.5323*(x>=.8)*(y<.8)  ...
    +.5323*(x<.8)*(y>=.8)  ...
    +1.5*(x>=.8)*(y>=.8);
u = 1.2060*(x<.8);
v = 1.2060*(y<.8);
E = p/(appdata.gamma-1)+rho*(u^2+v^2)/2;
qIC(4,1) = E;
qIC(3,1) = rho*v;
qIC(2,1) = rho*u;
qIC(1,1) = rho;


end

