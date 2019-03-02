function [ fq ] = problem_F(q,v1,v2)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
    %Shallow Water
    h = q(1);
    m = q(2);
    p = q(3);    
    alpha = 1/2;
    hinv = 1./h;
    fq =  [ m
    m.^2.*hinv+alpha.*h.^2
    m.*p.*hinv];
end



