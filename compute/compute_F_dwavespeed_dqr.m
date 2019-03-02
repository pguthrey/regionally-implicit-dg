function [ dspeeddq ] = compute_F_dwavespeed_dqr(ql,qr,quadpoint,data)
% Evaluates the data.numerical flux for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    ql
%           qr
% OUTdata.PUTS   Flux
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

qa = (ql+qr)/2;
[rhol] = problem_F_Jacobianspectralradius(ql,quadpoint,data);
[rhoa] = problem_F_Jacobianspectralradius(qa,quadpoint,data);
[rhor] = problem_F_Jacobianspectralradius(qr,quadpoint,data);
[~,imax] = max([rhol rhoa rhor]);

dspeeddq_tilde(3,:) = problem_F_drho_dq(qr,quadpoint,data);
dspeeddq_tilde(2,:) = problem_F_drho_dq(qa,quadpoint,data)/2;
dspeeddq_tilde(1,:) = problem_F_drho_dq(ql,quadpoint,data).*0;

dspeeddq = dspeeddq_tilde(imax,:);

end

