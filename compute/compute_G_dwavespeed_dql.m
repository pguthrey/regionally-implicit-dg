function [ dspeeddq ] = compute_G_dwavespeed_dql(ql,qr,quadpoint,data)
% Evaluates the data.numerical flux for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    ql
%           qr
% OUTdata.PUTS   Glux
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

qa = (ql+qr)/2;
[rhol] = problem_G_Jacobianspectralradius(ql,quadpoint,data);
[rhoa] = problem_G_Jacobianspectralradius(qa,quadpoint,data);
[rhor] = problem_G_Jacobianspectralradius(qr,quadpoint,data);
[~,imax] = max([rhol rhoa rhor]);

dspeeddq_tilde(3,:) = problem_G_drho_dq(qr,quadpoint,data).*0;
dspeeddq_tilde(2,:) = problem_G_drho_dq(qa,quadpoint,data)/2;
dspeeddq_tilde(1,:) = problem_G_drho_dq(ql,quadpoint,data);

dspeeddq = dspeeddq_tilde(imax,:);

end

