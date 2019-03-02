function [ Jacobian ] = problem_G_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
Rm = 1;
Rq = 1;
v1 = quadpoint(1);
v2 = quadpoint(2);
gamma = 1;%sqrt(Rm^2-v1.^2-v2.^2);

Jacobian = 0;%Rq*(appdata.E2-v1.*appdata.B3./gamma);

end

