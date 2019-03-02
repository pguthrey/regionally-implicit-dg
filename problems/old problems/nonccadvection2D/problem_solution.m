function [ exact ] = problem_solution(t,quadpoint,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
%load('data/data_parameters')
%load('data/data_mesh','v1centers','v2centers','deltav1','deltav2')

% Non-CC Advection
x = quadpoint(1);
y = quadpoint(2);
xprime = cos(t).*x-sin(t).*y;
yprime = sin(t).*x+cos(t).*y; 

exact = problem_IC([xprime yprime],data.appdata);
%exact= problem_IC(sqrt(x.^2+y.^2).*cos(atan2(y,x)-t),sqrt(x.^2+y.^2).*sin(atan2(y,x)-t));
%problem_IC(sqrt(x.^2+y.^2).*cos(atan2(y,x)-t),sqrt(x.^2+y.^2).*sin(atan2(y,x)-t));

end