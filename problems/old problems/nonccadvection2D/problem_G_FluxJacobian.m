function [ Jacobian ] = problem_G_FluxJacobian(q,quadpoint,appdata)
% Evaluates the flux Jacobian for the given problem
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    q
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

%CC Advection
%Jacobian = 1;
gamma = sqrt(1+quadpoint(1)^2+quadpoint(2)^2);
%NonCC Advection
Jacobian = -quadpoint(1)/gamma;%quadpoint(1);%*(max(abs(quadpoint(1:2))) < 1);

end

