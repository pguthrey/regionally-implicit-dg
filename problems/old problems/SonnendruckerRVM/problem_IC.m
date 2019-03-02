function [ qIC ] = problem_IC(point,appdata)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------re

x = point(1);
y = point(2);

%NonCC Advection
%{
qIC = exp(-200.*(sqrt(x.^2+y.^2)-.35)^2);
qIC = qIC + exp(-200.*((x-.35).^2+(y-0).^2));
qIC = qIC + exp(-200.*((x).^2+(y+.35).^2));
qIC = qIC + exp(-200.*((x).^2+(y-.35).^2));
%}
qIC = exp(-200.*((x+.35).^2+(y-0).^2));


end

