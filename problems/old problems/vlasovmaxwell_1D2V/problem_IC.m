function [ qIC ] = problem_IC(point,data)
% Evaluates the Initial Conditions for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y     : locations at which the ICs are to be evaluated
% OUTdata.PUTS   qIC     : ICs evaluated at the given locations
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

v1 = point(1);
v2 = point(2);
x = data.appdata.x;
fcoeffs = data.appdata.fcoeffs;
initial = @(v1,v2) DGeval_1D2V(fcoeffs,x,v1,v2,data);
qIC = initial(v1,v2);

end

