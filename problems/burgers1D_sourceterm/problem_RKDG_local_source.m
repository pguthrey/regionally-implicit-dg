function [rhs] = problem_RKDG_local_source(qstar,rhs,data,cellcenter)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INPUTS    DGsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTPUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

vectphi =  data.vectphi;
vectphi_Trans = data.vectphi_Trans;

nuv1 = data.nuv1;

Pd = data.Pd;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;

for k = 1:Pd
        v1quad = Dlist(k,1); 
        weight = Dquadwgts(k);
        wstar = vectphi(:,:,v1quad)*qstar;

        source = data.beta*wstar;
        
        fterm = (nuv1*weight)*vectphi_Trans(:,:,v1quad)*source;
        rhs = rhs - fterm;         
end
