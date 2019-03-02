function [ Jacobian ] = problem_F_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

u = q(2)/q(1);
N = u^2;
p = (appdata.gamma-1)*(q(3) - q(1)*N/2);
c = appdata.gamma*p/q(1);
h = (q(3)+p)/q(1);

Jacobian(3,:) = [-appdata.gamma*q(3)*u/q(1)+u^3*(appdata.gamma-1) ; ...
                    h - u^2*(appdata.gamma-1) ; ...
                    appdata.gamma*u ];
Jacobian(2,:) = [ (appdata.gamma-3)*u^2/2 ; ...
                    (3-appdata.gamma)*u ; ...
                    (appdata.gamma-1)];
Jacobian(1,:) = [0 1 0];

%{
x = quadpoint(1);
y = quadpoint(2);
check = (x>0)*(y>0)*(x<1)*(y<1);
Jacobian = Jacobian*check;
%}

end

