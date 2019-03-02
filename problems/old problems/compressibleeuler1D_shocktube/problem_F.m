function [ fq ] = problem_F(q,quadpoint,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

p = (appdata.gamma-1)*(q(3) - q(2)^2/(2*q(1)));
fq= [q(2); q(2)^2/q(1) + p  ; q(2)*(q(3)+p)/q(1) ];


end



