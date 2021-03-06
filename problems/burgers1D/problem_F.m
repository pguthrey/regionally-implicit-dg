function [ fq ] = problem_F(q,~,~)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------
error('whoops?')

p = 2/3*q(3) - q(2)^2/(3*q(1));
fq= [q(2); q(2)^2/q(1) + p  ; q(2)*(q(3)+p)/q(1) ];


end



