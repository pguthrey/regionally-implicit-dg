function [trunceast_temp,truncwest_temp, ...
          fluxeast_temp, fluxwest_temp, ...
          residual_cell,maxspeedF] ...
          = problem_precompute_eastwest_residual(qstar,qpast,qeast,qwest,data,~)
% Forms the %Jacobian of equations given by the DF data.method for a fixed cell
% written by data.Pierson Futhrey
% -------------------------------------------------
% INdata.PUTS    DFsubregion
%           prevcell
%           starcell
%           %Jac_truncdirs
% OUTdata.PUTS   %Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

residual_cell = (data.I_futr - data.Psi_dtau)*qstar-data.I_past*qpast;

trunceast_temp = 0; 
truncwest_temp = 0; 

fluxeast_temp = 0;
fluxwest_temp = 0;    

maxspeedF = 0;

nuv1 = data.nuv1;

vectpsi_east = data.vectpsi_east;
vectpsi_west = data.vectpsi_west;
vectpsi = data.vectpsi;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;
vectpsi_west_Trans = data.vectpsi_west_Trans;
vectpsi_east_Trans = data.vectpsi_east_Trans;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;
r_param = data.r_param;

for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    if space_dims >= 2
        v1quad = Dlist(k,2);
    else
        v1quad = 1;
    end
    if space_dims >= 3
        v2quad = Dlist(k,3);
    else
        v2quad = 1;
    end    
    weight = Dquadwgts(k);    
    wwest_p = vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest;
    weast_m = vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;
    wstar_p = vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    wstar_m = vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;

    Fwest_p = data.appdata.nuf*wwest_p;
    Feast_m = data.appdata.nuf*weast_m;

    Fstar_p = data.appdata.nuf*wstar_p;    
    Fstar_m = data.appdata.nuf*wstar_m;
    
    speedeast = abs(data.appdata.nuf);
    speedwest = abs(data.appdata.nuf);
        
    %east Boundary terms
    psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    trunceast_temp = trunceast_temp + (psik*Fstar_p)*(nuv1*weight);
    if r_param > 0
        Flux = (Fstar_p + Feast_m - speedeast*(weast_m-wstar_p))/2;%  compute_F_numerical_flux(wstar,weast,quadpoint,appdata);
        fluxeast_temp = fluxeast_temp + (psik*Flux)*(nuv1*weight);  
    end
    
    %west boundary terms

    psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    truncwest_temp = truncwest_temp - (psik*Fstar_m)*(nuv1*weight);
    if r_param > 0
        Flux =  (Fstar_m + Fwest_p - speedwest*(wstar_m-wwest_p))/2;%compute_F_numerical_flux(wwest,wstar,quadpoint,appdata);
        fluxwest_temp = fluxwest_temp - (psik*Flux)*(nuv1*weight);
    end
    
    maxspeedF_quad = max(speedeast,speedwest);
    maxspeedF = max(maxspeedF,maxspeedF_quad);
    
end

for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        if space_dims >= 2
            v2quad = Dp1list(k,3);
        else
            v2quad = 1;
        end
        if space_dims >= 3
            v3quad = Dp1list(k,4);
        else
            v3quad = 1;
        end
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        Fstar = data.appdata.nuf*wstar;
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        gterm = (nuv1*weight)*psikdxii*Fstar;
        residual_cell = residual_cell - gterm; 
        maxspeedF = max(maxspeedF,abs(data.appdata.nuf);
end


end