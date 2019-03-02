function [ Jacobian ] = compute_F_dflux_dql(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dFdql = problem_F_FluxJacobian(ql,quadpoint,data);
[speed] = compute_F_wavespeed(ql,qr,quadpoint,data);
[dspeeddql] = 0;%compute_F_dwavespeed_dql(ql,qr,quadpoint,data);

Jacobian = (dFdql+speed.*(eye(length(ql))) - (qr-ql)*dspeeddql)./2;

end

