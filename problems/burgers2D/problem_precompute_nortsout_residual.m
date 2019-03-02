function [truncnort_temp,truncsout_temp, ...
          fluxnort_temp, fluxsout_temp, ...
          residual_cell,maxspeedG] ...
        = predictor_precompute_nortsout_residual(qstar,qnort,qsout,residual_cell,data,cellcenter)
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

speednort2 = 0;
speedsout2 = 0;

v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);

for k = 1:data.Pd
    tquad = data.Dlist(k,1); 
    %tloc = data.Dquadlocs(k,1);
    if data.space_dims >= 2
        v1quad = data.Dlist(k,2);
        v1loc = data.Dquadlocs(k,2);
    else
        v1quad = 1;
        v1loc = inf;
    end
    if data.space_dims >= 3
        v2quad = data.Dlist(k,3);
        v2loc = data.Dquadlocs(k,3);
    else
        v2quad = 1;
        v2loc = inf;
    end    
    weight = data.Dquadwgts(k);    
    wsout = data.vectpsi_nort(:,:,tquad,v1quad,v2quad)*qsout;
    wnort = data.vectpsi_sout(:,:,tquad,v1quad,v2quad)*qnort;

    %North Boundary terms
    v1 = v1center+data.deltav1/2*v1loc;        
    v2 = v2center+data.deltav2/2;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];
    wstar = data.vectpsi_nort(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    [ Flux ] = problem_G(wstar,quadpoint,data.appdata);
    [speednort1] = problem_G_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    truncnort_temp = truncnort_temp + data.nuv2*(psik*Flux)*weight;
    if data.r_param > 0
        [ Flux,speednort2] = compute_G_numerical_flux(wstar,wnort,quadpoint,data.appdata);
        fluxnort_temp = fluxnort_temp + data.nuv2*(psik*Flux)*weight;  
    end
    
    %sout boundary terms
    v1 = v1center-data.deltav1/2*v1loc;        
    v2 = v2center+data.deltav2/2;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];

    wstar = data.vectpsi_sout(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    [Flux] = problem_G(wstar,quadpoint,data.appdata);
    [speedsout1] = problem_G_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    truncsout_temp = truncsout_temp - data.nuv2*(psik*Flux)*weight;
    if data.r_param > 0
        [ Flux,speedsout2 ] = compute_G_numerical_flux(wsout,wstar,quadpoint,data.appdata);
        fluxsout_temp = fluxsout_temp - data.nuv2*(psik*Flux)*weight;
    end
    
    maxspeedG = max([speednort1 speednort2 speedsout1 speedsout2]);
end

v2center = 0;
v3center = 0;

for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        %tloc = data.Dp1quadlocs(k,1);
        v1quad = data.Dp1list(k,2); 
        v1loc = data.Dp1quadlocs(k,2);
        v1center = cellcenter(1);
        if data.space_dims >= 2
            v2quad = data.Dp1list(k,3);
            v2loc = data.Dp1quadlocs(k,3);
            v2center = cellcenter(2);
        else
            v2quad = 1;
            v2loc = inf;
        end
        if data.space_dims >= 3
            v3quad = data.Dp1list(k,4);
            v3loc = data.Dp1quadlocs(k,4);
            v3center = cellcenter(3);
        else
            v3quad = 1;
            v3loc = inf;
        end
        weight = data.Dp1quadwgts(k);
        %Derivative terms
        v1 = v1center+data.deltav1/2*v1loc;
        v2 = v2center+data.deltav2/2*v2loc;
        v3 = v3center+data.deltav3/2*v3loc;
        quadpoint = [v1 v2 v3];
        wstar = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        psikdeta = data.vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        gterm = (data.nuv2*weight)*psikdeta*problem_G(wstar,quadpoint,data.appdata);
        residual_cell = residual_cell - gterm; 
end
end