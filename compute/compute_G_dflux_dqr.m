function [ Jacobian ] = compute_G_dflux_dqr(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dGdqr = problem_G_FluxJacobian(qr,quadpoint,data);
[speed] = compute_G_wavespeed(ql,qr,quadpoint,data);
[dspeeddqr] = compute_G_dwavespeed_dqr(ql,qr,quadpoint,data);

Jacobian = (dGdqr-speed.*(eye(length(qr))) - (qr-ql)*dspeeddqr)./2;

end

