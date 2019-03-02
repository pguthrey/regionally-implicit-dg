function [ Flux ,speed ] = compute_F_numerical_flux(ql,qr,quadpoint,data)
% Evaluates the data.numerical flux for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    ql
%           qr
% OUTdata.PUTS   Flux
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

%Local Lax-Friedrichs
Fl = problem_F(ql,quadpoint,data);
Fr = problem_F(qr,quadpoint,data);

[speed] = compute_F_wavespeed(ql,qr,quadpoint,data);

Flux = (Fl+Fr-speed.*(qr-ql))./2;

end

