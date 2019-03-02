function [residual] ...
          = problem_res_init(qstar,qpast,data)
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

vectpsi_dxii_Trans  = data.vectpsi_dxii_Trans;
vectpsi_deta_Trans  = data.vectpsi_deta_Trans;
vectpsi = data.vectpsi;

v2quad = 1;
v3quad = 1;
Dp1list = data.Dp1list;
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1quadwgts = data.Dp1quadwgts;
nuv1 = data.nuv1;
nuv2 = data.nuv2;

%%{
% RES ONLY 
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
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        Fstar = wstar.^2/2;
        fterm = (nuv1*weight)*psikdxii*Fstar;
        residual = residual  - fterm;
        
        psikdeta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        Gstar = wstar.^2/2;
        gterm = (nuv2*weight)*psikdeta*Gstar;
        residual = residual  - gterm;       
end


% JAC ONLY
%{
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
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);        
        DFstardq = wstar;
        Jacobian = Jacobian - (nuv1*weight)*psikdxii*DFstardq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);

        psikdeta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);        
        DGstardq = wstar;
        Jacobian = Jacobian - (nuv2*weight)*psikdeta*DGstardq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
end
%}








end