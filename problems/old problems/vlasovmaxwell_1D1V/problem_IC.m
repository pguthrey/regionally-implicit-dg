function [ qIC ] = problem_IC(point,data)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

v1 = point(1);
x = data.appdata.x;
fcoeffs = data.appdata.fcoeffs;
initial = @(v1) DGeval_1D1V(fcoeffs,x,v1,data);
qIC = initial(v1);

end

