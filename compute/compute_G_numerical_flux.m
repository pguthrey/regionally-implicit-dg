function [ Glux ,speed ] = compute_G_numerical_flux(ql,qr,quadpoint,data)
% Evaluates the data.numerical flux for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    ql
%           qr
% OUTdata.PUTS   Glux
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

%Local Lax-Griedrichs
qa = (ql+qr)/2;
Gl = problem_G(ql,quadpoint,data);
Gr = problem_G(qr,quadpoint,data);

[speed] = compute_G_wavespeed(ql,qr,quadpoint,data);

Glux = (Gl+Gr-speed.*(qr-ql))./2;

end

