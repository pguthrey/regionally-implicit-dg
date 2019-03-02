function [ fq ] = problem_F(q,quadpoint,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

    %CC Advection
    %fq = q;
    
    %Non-CC advection
    fq= v2*q;%-quadpoint(2)*q;%*(max(abs(quadpoint(1:2))) < 1);
end



