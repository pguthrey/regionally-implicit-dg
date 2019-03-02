function [ Jacobian ] = problem_G_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

Jacobian =  q;

%{
x = quadpoint(1);
y = quadpoint(2);

check = (x>0)*(y>0)*(x<1)*(y<1);
Jacobian = Jacobian*check;
%}

end

