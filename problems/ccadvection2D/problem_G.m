function [ gq ] = problem_G(q,quadpoint,appdata)
% Evaluates g(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   g evaluated for conserved variables g
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

%CC Advection
%    gq = q;
    
%NON-CC Advection
gq = appdata.nug*q;

end




