function [ qIC ] = problem_IC(quadpoint,appdata)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------re

%NonCC Advection
x = quadpoint(1);
qIC = .5+sin(3*pi*x)/2;
qIC = (x-.3).^3.*(x+.3).^3;
qIC = exp(-200.*(x.^2));
qIC = .5+sin(3*pi*x)/2;
qIC = exp(-300.*(x.^2));
qIC = exp(-200.*(x.^2));
qIC = cos(10*pi*x);
qIC = sin(pi*x);
qIC = exp(-300.*(x.^2));
qIC = exp(-300.*(x.^2));
qIC = sin(pi*x);
qIC = cos(4*pi*x);
qIC = exp(-400.*(x).^4).*x^2.*20.*cos(x/10)+sin(x).*exp(-200*x.^2);
omega = 4;
qIC = cos(omega*pi*x);
%qIC = .3;
end

