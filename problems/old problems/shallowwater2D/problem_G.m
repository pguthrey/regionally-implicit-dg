function [ gq ] = problem_G(q,v1,v2)
% Evaluates g(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   g evaluated for conserved variables g
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------


%Shallow Water
h = q(1);
m = q(2);
p = q(3);

alpha = 1/2;

hinv = 1./h;
gq =  [ p
        m.*p.*hinv
        p.^2.*hinv+alpha.*h.^2];
end




