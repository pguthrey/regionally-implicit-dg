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

qIC = .3+0.*exp(-200.*((x+.25).^2));
qIC = .5+0.*x;
qIC = exp(-200.*((x+.25).^2));
qIC = .4*(abs(x) < .2 );
end

