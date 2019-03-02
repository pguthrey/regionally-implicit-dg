function [ fq ] = problem_F(q,quadpoint,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

rho = q(1);
u = q(2)/q(1);
v = q(3)/q(1);
E = q(4);
speedsq = u^2+v^2;

p = (appdata.gamma-1)*(E - speedsq*rho/2);
fq= [ rho*u; rho*u^2 + p ; rho*u*v ; u*(E+p) ];

%{
p = (appdata.gamma-1)*(q(4) - (q(2)^2+q(3)^2)/(2*q(1)));
fq= [q(2); q(2)^2/q(1) + p ; q(2)*q(3)/q(1) ; q(2)*(q(4)+p)/q(1) ];
%}

end



