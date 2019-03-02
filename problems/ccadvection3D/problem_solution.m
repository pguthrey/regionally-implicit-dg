function [ exact ] = problem_solution(time,point,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2,v3) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

t = time;
v1 = point(1);
v2 = point(2);
v3 = point(3);

%CC Advection
arg1 = v1-data.appdata.nuf*t;
arg2 = v2-data.appdata.nug*t;
arg3 = v3-data.appdata.nuh*t;

perv1 = floor((arg1+data.v1_lb)/(data.v1_ub-data.v1_lb))+1;
perv2 = floor((arg2+data.v2_lb)/(data.v2_ub-data.v2_lb))+1;
perv3 = floor((arg3+data.v3_lb)/(data.v3_ub-data.v3_lb))+1;

v1tilde = arg1-(data.v1_ub-data.v1_lb)*perv1;
v2tilde = arg2-(data.v2_ub-data.v2_lb)*perv2;
v3tilde = arg3-(data.v3_ub-data.v3_lb)*perv3;

exact = problem_IC([v1tilde,v2tilde,v3tilde]);

end