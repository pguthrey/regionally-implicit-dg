function [ Jacobian ] = compute_G_dflux_dql(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dGdql = problem_G_FluxJacobian(ql,quadpoint,data);
[speed] = compute_G_wavespeed(ql,qr,quadpoint,data);
[dspeeddql] = compute_G_dwavespeed_dql(ql,qr,quadpoint,data);

Jacobian = (dGdql+speed.*(eye(length(ql))) - (qr-ql)*dspeeddql)./2;

end

