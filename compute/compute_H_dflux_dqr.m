function [ Jacobian ] = compute_H_dflux_dqr(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dHdqr = problem_H_FluxJacobian(qr,quadpoint,data);
[speed] = compute_H_wavespeed(ql,qr,quadpoint,data);
[dspeeddqr] = compute_H_dwavespeed_dqr(ql,qr,quadpoint,data);

Jacobian = (dHdqr-speed.*(eye(length(qr))) - (qr-ql)*dspeeddqr)./2;

end

