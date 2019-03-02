function [ gq ] = problem_G(q,quadpoint,appdata)
% Evaluates g(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   g evaluated for conserved variables g
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

%CC Advection
%    gq = q;
    
%NON-CC Advection
gq = quadpoint(1)*q;%*(max(abs(quadpoint(1:2))) < 1);

end




