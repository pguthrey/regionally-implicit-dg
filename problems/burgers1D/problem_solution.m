function [ exact ] = problem_solution(t,quadpoint,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2,v3) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

u0 = problem_IC(quadpoint,data.appdata);
x = quadpoint(1);


options = optimoptions('fsolve','Display','none');
options.FunctionTolerance = 1e-14;
options.OptimalityTolerance = 1e-14;
options.StepTolerance = 1e-14;


exact = fsolve(@(u)residual(u,x,t),u0,options);


end


function [res] = residual(u,x,t)

    res = u - problem_IC(x-u*t,NaN);

end