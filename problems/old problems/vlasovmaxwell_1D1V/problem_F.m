function [ fq ] = problem_F(q,quadpoint,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
Rq = -1;
v1 = quadpoint(1);
gamma = 1;%sqrt(Rm^2-v1.^2-v2.^2);
check = (v1 < appdata.v1_ub)*(v1 > appdata.v1_lb);
fq = Rq*(appdata.E1)*q*check;
end



