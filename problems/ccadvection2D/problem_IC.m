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
qIC = exp(-200.*((x).^2+(y).^2));
%qIC = exp(-200.*((y+.5).^2));
%qIC = sin(pi*x).*sin(pi*y);%exp(-200.*((y).^2));
%omega = 16;
%qIC = sin(omega*pi*x).*sin(omega*pi*y);%exp(-200.*((y).^2));
%qIC = .5;
%qIC = exp(-200.*((x).^2+(y).^2));
end

