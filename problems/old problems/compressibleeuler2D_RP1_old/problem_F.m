function [ fq ] = problem_F(q,quadpoint,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

p = (appdata.gamma-1)*(q(4) - (q(2)^2+q(3)^2)/(2*q(1)));
fq= [q(2); q(2)^2/q(1) + p ; q(2)*q(3)/q(1) ; q(2)*(q(4)+p)/q(1) ];

%{
x = quadpoint(1);
y = quadpoint(2);

check = (x>0)*(y>0)*(x<1)*(y<1);
fq = fq*check;
%}

end



