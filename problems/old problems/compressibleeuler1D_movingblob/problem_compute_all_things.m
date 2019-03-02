function [F,Jacobian,spd,dspeeddq] = problem_compute_all_things(q,~,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------
m = q(2);
orho = 1/q(1);
u = q(2)*orho;
u2 = u*u;
u3 = u2*u;
gm1 = appdata.gamma-1;
gm3 = gm1-2;
p = gm1*(q(3) - q(1)*u2*0.5);
h = u*(q(3)+p);
c = sqrt(appdata.gamma*p*orho);
%{
Jacobian = zeros(3); 

Jacobian(3,:) = [ -appdata.gamma*u*q(3)*orho + gm1*u3 ; ...
                    appdata.gamma*q(3)*orho - gm1*1.5*u2 ; ...
                    appdata.gamma*u ];
Jacobian(2,:) = [ gm3*u2*0.5 ; ...
                    -gm3*u ; ...
                    gm1];
Jacobian(1,:) = [0 1 0];
%}
Jacobian = [[0 1 0];
            [ gm3*u2*0.5 , ...
                    -gm3*u , ...
                    gm1];
            [ -appdata.gamma*u*q(3)*orho + gm1*u3 , ...
                    appdata.gamma*q(3)*orho - gm1*1.5*u2 , ...
                    appdata.gamma*u ];];
                
F = [m ; m*u + p  ; h ];

spd = abs(u)+c;

%dcdq = appdata.gamma*gm1/(2*c*p^2)*[E -m -rho];
dspeeddq = 0;%[abs(u) sign(u) 0]*orho + dcdq;

end



