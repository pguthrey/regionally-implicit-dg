function [ Jacobian ] = problem_F_FluxJacobian(q,~,~)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

u = q(2)/q(1);
u2 = u*u;
u3 = u2*u;
g = 5/2;
gm1 = 2/3;
gm3 = gm1-2;
p = gm1*(q(3) - q(1)*u2/2);
h = q(2)*(q(3)+p)/q(1);

Jacobian = zeros(3); 
Jacobian(3,:) = [-g*q(3)*u/q(1)+u^3*gm1 ; ...
                    h - u^2*gm1 ; ...
                    g*u ];
Jacobian(2,:) = [ gm3*u2/2 ; ...
                    -gm3*u ; ...
                    gm1];
Jacobian(1,:) = [0 1 0];

u = q(2)/q(1);
N = u^2;
p = gm1*(q(3) - q(1)*N/2);
c = g*p/q(1);
h = (q(3)+p)/q(1);

Jacobian(3,:) = [-g*q(3)*u/q(1)+u^3*gm1 ; ...
                    h - u^2*gm1 ; ...
                    g*u ];
Jacobian(2,:) = [ gm3*u^2/2 ; ...
                    -gm3*u ; ...
                    gm1];
Jacobian(1,:) = [0 1 0];

%{
x = quadpoint(1);
y = quadpoint(2);
check = (x>0)*(y>0)*(x<1)*(y<1);
Jacobian = Jacobian*check;
%}

end

