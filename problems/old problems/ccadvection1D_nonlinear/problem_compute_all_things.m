function [F,Jacobian,spd,dspeeddq] = problem_compute_all_things(q,~,appdata)
% Evaluates f(q) for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   f evaluated for conserved variables f
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

F = appdata.nuf*q;
Jacobian = appdata.nuf;
spd = abs(appdata.nuf);
dspeeddq = 0;

end



