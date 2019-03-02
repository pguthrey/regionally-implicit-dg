function [ exact ] = problem_solution(t,quadpoint,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2,v3) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

v1 = quadpoint(1);
%CC Advection
arg1 = v1-data.appdata.nuf*t;

perv1 = floor((arg1+data.v1_lb)/(data.v1_ub-data.v1_lb))+1;

v1tilde = arg1-(data.v1_ub-data.v1_lb)*perv1;

exact = problem_IC(v1tilde);

end