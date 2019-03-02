function [ gq ] = problem_G(q,quadpoint,appdata)
% Evaluates g(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   g evaluated for conserved variables g
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

p = (appdata.gamma-1)*(q(4) - (q(2)^2+q(3)^2)/(2*q(1)));
gq= [q(3); q(3)*q(2)/q(1)  ; q(3)^2/q(1) + p ; q(3)*(q(4)+p)/q(1) ];

%{
x = quadpoint(1);
y = quadpoint(2);
%}
%{
check = (x>0)*(y>0)*(x<1)*(y<1);
gq = gq*check;
%}

end




