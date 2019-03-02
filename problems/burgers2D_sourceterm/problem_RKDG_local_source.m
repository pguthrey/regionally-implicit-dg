function [residual] = problem_RKDG_local_source(qstar,residual,data,cellcenter)
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

v2quad = 1;
v3quad = 1;
rhs_cell = 0;

Pd = data.Pd;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;

space_dims = data.space_dims;

for k = 1:Pd
    v1quad = Dlist(k,1);
    if space_dims >= 2
        v2quad = Dlist(k,2);
        if space_dims >= 3
            v3quad = Dlist(k,3);
        end
    end
    weight = Dquadwgts(k);
    wstar = vectphi(:,:,v1quad,v2quad,v3quad)*qstar;
    phi = vectphi(:,:,v1quad,v2quad,v3quad)';
    u = wstar(1);            
    Source(1,1) = data.beta*u ;
    rhs_cell = rhs_cell - phi*Source*(weight);
end


residual = residual + rhs_cell;
