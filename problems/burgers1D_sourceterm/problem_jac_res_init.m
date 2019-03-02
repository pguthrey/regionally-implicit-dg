function [Jacobian,residual] ...
          = problem_jac_res_init(qstar,qpast,data)
% Forms the %Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           %Jac_truncdirs
% OUTdata.PUTS   %Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

residual = (data.I_futr - data.Psi_dtau)*qstar-data.I_past*qpast;
Jacobian = (data.I_futr - data.Psi_dtau);

vectpsi_dxii_Trans  = data.vectpsi_dxii_Trans;
vectpsi = data.vectpsi;

v2quad = 1;
v3quad = 1;
Dp1list = data.Dp1list;
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1quadwgts = data.Dp1quadwgts;
nuv1 = data.nuv1;

for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            if space_dims >= 3
                v3quad = Dp1list(k,4);
            end
        end
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        Fstar = wstar.^2/2;
        DFstar_dq = wstar;
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        fterm = (nuv1*weight)*psikdxii*Fstar;
        residual = residual  - fterm;
        Jacobian = Jacobian - (nuv1*weight)*psikdxii*DFstar_dq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);

end



end