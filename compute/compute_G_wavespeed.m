function [ speed ] = compute_G_wavespeed(ql,qr,quadpoint,data)
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

speed = max([rhol rhoa rhor]);

end

