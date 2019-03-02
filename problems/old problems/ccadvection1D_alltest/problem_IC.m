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
qIC = (x>-0.35)*(x<-0.15);
qIC = qIC + (x>-0.85)*(x<-0.65)*cos((x+.75)/.1*pi/2);
qIC = qIC + (x>.15)*(x<.35)*(1-abs(x-.25)/.1);
qIC = qIC + (x>.65)*(x<.85)*exp(-1e3.*(x-.75)^2);

end

