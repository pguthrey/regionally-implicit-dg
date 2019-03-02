function [ Jacobian ] = problem_G_FluxJacobian(q,v1,v2)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

%Shallow water
h = q(1);
m = q(2);
p = q(3);
alpha = 1/2;
u = m./h;
v = p./h;
Jacobian = [[0,0,1];
            [-u.*v,v,u];
            [-v.^2+2.*alpha.*h,0,2.*v]];

end

