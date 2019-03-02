function [ Jacobian ] = problem_F_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
Rq = -1;
v1 = quadpoint(1);
gamma = 1;%sqrt(Rm^2-v1.^2-v2.^2);
check = (v1 < appdata.v1_ub)*(v1 > appdata.v1_lb);
Jacobian = Rq*(appdata.E1)*check;
end

