function [ Jacobian ] = compute_H_dflux_dql(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dHdql = problem_H_FluxJacobian(ql,quadpoint,data);
[speed] = compute_H_wavespeed(ql,qr,quadpoint,data);
[dspeeddql] = compute_H_dwavespeed_dql(ql,qr,quadpoint,data);

Jacobian = (dHdql+speed.*(eye(length(ql))) - (qr-ql)*dspeeddql)./2;

end

