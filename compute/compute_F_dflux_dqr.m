function [ Jacobian ] = compute_F_dflux_dqr(ql,qr,quadpoint,data)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

dFdqr = problem_F_FluxJacobian(qr,quadpoint,data);
[speed] = compute_F_wavespeed(ql,qr,quadpoint,data);
[dspeeddqr] = 0;%compute_F_dwavespeed_dqr(ql,qr,quadpoint,data);

Jacobian = (dFdqr-speed.*(eye(length(qr))) - (qr-ql)*dspeeddqr)./2;

end

