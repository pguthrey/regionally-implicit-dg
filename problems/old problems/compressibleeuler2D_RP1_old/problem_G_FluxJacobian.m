function [ Jacobian ] = problem_G_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

u = q(2)/q(1);
v = q(3)/q(1);
N = u^2+v^2;
p = (appdata.gamma-1)*(q(4) - q(1)*N/2);
h = (q(4)+p)/q(1);

Jacobian(4,:) = [-appdata.gamma*q(4)*v/q(1)+v*N*(appdata.gamma-1), ...
                    - u*v*(appdata.gamma-1) , ...
                    h  - v^2*(appdata.gamma-1) , ...
                    appdata.gamma*v ];
Jacobian(3,:) = [-v^2+(appdata.gamma-1)*N/2 , ...
                    -(appdata.gamma-1)*u , ...
                    2*v-(appdata.gamma-1)*v , ... 
                    (appdata.gamma-1)];
Jacobian(2,:) = [-v*u , ...
                    v  , ...
                    u ,  ...
                    0];
Jacobian(1,:) = [0 0 1 0];

%{
x = quadpoint(1);
y = quadpoint(2);

check = (x>0)*(y>0)*(x<1)*(y<1);
Jacobian = Jacobian*check;
%}

end

