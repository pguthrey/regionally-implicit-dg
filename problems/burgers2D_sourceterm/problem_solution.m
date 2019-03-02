function [ exact ] = problem_solution(t,quadpoint,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2,v3) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------


%%{
u0 = problem_IC(quadpoint,data.appdata);
x = quadpoint(1);
y = quadpoint(2);


%gamma = 


options = optimoptions('fsolve','Display','none');
options.FunctionTolerance =  4*eps;
options.OptimalityTolerance =  4*eps;
options.StepTolerance =  4*eps;

exact = fsolve(@(u)residual(u,x,y,t),u0,options);

%}

%{
v1 = quadpoint(1);
v2 = quadpoint(2);

%CC Advection
arg1 = v1-data.appdata.xspd*t;
arg2 = v2-data.appdata.yspd*t;

perv1 = floor((arg1+data.v1_lb)/(data.v1_ub-data.v1_lb))+1;
perv2 = floor((arg2+data.v2_lb)/(data.v2_ub-data.v2_lb))+1;

v1tilde = arg1-(data.v1_ub-data.v1_lb)*perv1;
v2tilde = arg2-(data.v2_ub-data.v2_lb)*perv2;

exact = problem_IC([v1tilde,v2tilde]);
%}

end


function [res] = residual(u,x,y,t)

    xtilde = x-u*t;
    ytilde = y-u*t;
    res = u - problem_IC([xtilde,ytilde],NaN);

end