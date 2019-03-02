function [ Jacobian ] = problem_F_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

%CC Advection
%Jacobian = 1;

%NonCC Advection
Jacobian = appdata.nuf;


end

