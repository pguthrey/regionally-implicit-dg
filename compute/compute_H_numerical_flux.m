function [ Flux ,speed ] = compute_H_numerical_flux(ql,qr,quadpoint,data)
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
Hl = problem_H(ql,quadpoint,data);
Hr = problem_H(qr,quadpoint,data);

[speed] = compute_H_wavespeed(ql,qr,quadpoint,data);

Flux = (Hl+Hr-speed.*(qr-ql))./2;

end

