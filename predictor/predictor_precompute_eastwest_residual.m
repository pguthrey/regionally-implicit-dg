function [trunceast_temp,truncwest_temp, ...
          fluxeast_temp, fluxwest_temp, ...
          residual_cell,maxspeedF] ...
          = predictor_precompute_eastwest_residual(qstar,qpast,qeast,qwest,data,cellcenter)
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

trunceast_temp = 0; 
truncwest_temp = 0; 

fluxeast_temp = 0;
fluxwest_temp = 0;    

residual_cell = (data.I_futr - data.Psi_dtau)*qstar-data.I_past*qpast;

maxspeedF = 0;

v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);

speedeast2 = 0;
speedwest2 = 0;
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
    wwest = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest;
    weast = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;

    %East Boundary terms
    v1 = v1center+data.deltav1/2;        
    v2 = v2center+data.deltav2/2*v1loc;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];
    wstar = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    [ Flux ] = problem_F(wstar,quadpoint,data.appdata);
    [speedeast1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    trunceast_temp = trunceast_temp + data.nuv1*(psik*Flux)*weight;
    if data.r_param > 0
        [ Flux,speedeast2] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
        fluxeast_temp = fluxeast_temp + data.nuv1*(psik*Flux)*weight;  
    end
    
    %west boundary terms
    v1 = v1center-data.deltav1/2;
    v2 = v2center+data.deltav2/2*v1loc;
    v3 = v3center+data.deltav3/2*v2loc;
    quadpoint = [v1 v2 v3];

    wstar = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    [Flux] = problem_F(wstar,quadpoint,data.appdata);
    [speedwest1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    truncwest_temp = truncwest_temp - data.nuv1*(psik*Flux)*weight;
    if data.r_param > 0
        [ Flux,speedwest2 ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
        fluxwest_temp = fluxwest_temp - data.nuv1*(psik*Flux)*weight;
    end
    maxspeedF = max([speedeast1 speedeast2 speedwest1 speedwest2 maxspeedF]);
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
        psikdxii = data.vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        fterm = (data.nuv1*weight)*psikdxii*problem_F(wstar,quadpoint,data.appdata);
        residual_cell = residual_cell  - fterm;
        [speedjac] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        maxspeedF = max([speedjac maxspeedF]);
end



end