function [truncnort_temp,truncsout_temp, ...
          fluxnort_temp, fluxsout_temp, ...
          residual_cell,maxspeedG] ...
        = problem_precompute_nortsout_residual(qstar,qnort,qsout,residual_cell,data,~)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           %truncdirs
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

truncnort_temp = 0; 
truncsout_temp = 0; 

fluxnort_temp = 0;
fluxsout_temp = 0;    

maxspeedG = 0;

nuv2 = data.nuv2;

vectpsi_nort = data.vectpsi_nort;
vectpsi_sout = data.vectpsi_sout;
vectpsi = data.vectpsi;
vectpsi_deta_Trans = data.vectpsi_deta_Trans;
vectpsi_sout_Trans = data.vectpsi_sout_Trans;
vectpsi_nort_Trans = data.vectpsi_nort_Trans;

Pd = data.Pd;
Pdp1 = data.Pdp1;
r_param = data.r_param;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

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
    wsout_p = vectpsi_nort(:,:,tquad,v1quad,v2quad)*qsout;
    wnort_m = vectpsi_sout(:,:,tquad,v1quad,v2quad)*qnort;
    wstar_p = vectpsi_nort(:,:,tquad,v1quad,v2quad)*qstar;
    wstar_m = vectpsi_sout(:,:,tquad,v1quad,v2quad)*qstar;

    Fsout_p = wsout_p^2/2;
    Fnort_m = wnort_m^2/2;

    Fstar_p = wstar_p^2/2;    
    Fstar_m = wstar_m^2/2;
    
    speednort = max(abs(wstar_p),abs(wnort_m));
    speedsout = max(abs(wstar_m),abs(wsout_p));
        
    %North Boundary terms
    psik = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    truncnort_temp = truncnort_temp + (psik*Fstar_p)*(nuv2*weight);
    if r_param > 0
        Flux = (Fstar_p + Fnort_m - speednort*(wnort_m-wstar_p))/2;%  compute_G_numerical_flux(wstar,wnort,quadpoint,appdata);
        fluxnort_temp = fluxnort_temp + (psik*Flux)*(nuv2*weight);  
    end
    
    %sout boundary terms

    psik = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    truncsout_temp = truncsout_temp - (psik*Fstar_m)*(nuv2*weight);
    if r_param > 0
        Flux =  (Fstar_m + Fsout_p - speedsout*(wstar_m-wsout_p))/2;%compute_G_numerical_flux(wsout,wstar,quadpoint,appdata);
        fluxsout_temp = fluxsout_temp - (psik*Flux)*(nuv2*weight);
    end
    
    maxspeedG_quad = max(speednort,speedsout);
    maxspeedG = max(maxspeedG,maxspeedG_quad);
    
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
        Fstar = wstar^2/2;
        psikdeta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        gterm = (nuv2*weight)*psikdeta*Fstar;
        residual_cell = residual_cell - gterm; 
        maxspeedG = max(maxspeedG,abs(wstar));

end
end