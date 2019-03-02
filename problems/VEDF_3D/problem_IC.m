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
z = point(3);

n = exp(-200.0*(x.^2+y.^2+z.^2))+.1;

qIC(:,1) = [n ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];

end
