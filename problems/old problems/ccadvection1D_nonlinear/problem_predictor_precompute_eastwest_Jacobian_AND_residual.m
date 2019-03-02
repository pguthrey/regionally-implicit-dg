function [Jac_east_cell_temp,Jac_west_cell_temp, ...
        Jac_east_other_temp,Jac_west_other_temp, ...
        Jac_trunc_east_temp,Jac_trunc_west_temp, ...
        Jac_cell, ...
        trunceast_temp,truncwest_temp, ...
        fluxeast_temp, fluxwest_temp, ...
        residual_cell,maxspeedF] ...
        = problem_predictor_precompute_eastwest_Jacobian_AND_residual(qstar,qeast,qwest,qpast,data,cellcenter)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

Jac_cell = data.I_futr - data.Psi_dtau;

Jac_east_cell_temp = zeros(data.thetaT,data.thetaT);
Jac_west_cell_temp = zeros(data.thetaT,data.thetaT);    

Jac_east_other_temp = zeros(data.thetaT,data.thetaT);
Jac_west_other_temp = zeros(data.thetaT,data.thetaT);  

Jac_trunc_east_temp = zeros(data.thetaT,data.thetaT);   
Jac_trunc_west_temp = zeros(data.thetaT,data.thetaT); 

trunceast_temp = 0; 
truncwest_temp = 0; 

fluxeast_temp = 0;
fluxwest_temp = 0;    

residual_cell = (data.I_futr - data.Psi_dtau)*qstar-data.I_past*qpast;

maxspeedF = 0;


speedeast2 = 0;
speedwest2 = 0;

%vectpsi_all = data.vectpsi_all;
vectpsi_east = data.vectpsi_east;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west = data.vectpsi_west;
vectpsi_west_Trans = data.vectpsi_west_Trans;
vectpsi = data.vectpsi;




nuv1 = data.nuv1;

vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;

nuf = data.appdata.nuf;


I_mat = eye(data.Neqns);
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
    wwest = vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest;
    weast = vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;

    %wwest = vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest;
    %weast = vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;
    
    %wwest2 = vectpsi_all(:,:,1+tquad,end)*qwest;
    %weast2 = vectpsi_all(:,:,1+tquad,1)*qeast;
    
    
    %err = wwest-wwest2
    %keyoard
    
    %East Boundary terms
    wstar = vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
       
    %[F_star,JF_star,speedF_star,dspeeddq_star] = problem_compute_all_things(wstar,quadpoint,data.appdata);
    F_star = nuf*wstar; 
    JF_star = nuf; 
    speedF_star = abs(nuf); 
    dspeeddq_star = 0;
        %F_star = problem_F(wstar,quadpoint,data.appdata);
    
    trunceast_temp = trunceast_temp + nuv1*(psik*F_star)*weight;
    if data.r_param > 0
        %[F_east,JF_east,speedF_east,dspeeddq_east] = problem_compute_all_things(weast,quadpoint,data.appdata);
        F_east = nuf*weast; 
        JF_east = nuf; 
        speedF_east = abs(nuf); 
        dspeeddq_east = 0;
        %[speedF_avg] = problem_F_Jacobianspectralradius((weast+wstar)*0.5,quadpoint,data.appdata);
        speedF_avg = abs(nuf);
        speedeast = max([speedF_star speedF_east speedF_avg]);
        Flux = (F_star+F_east-speedeast.*(weast-wstar))*0.5;
        %[ ~ , speedeast2] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
        %???
        fluxeast_temp = fluxeast_temp + nuv1*(psik*Flux)*weight;  
    end
    
    maxspeedF = max(speedeast2,maxspeedF);
    
    psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
    %JF_star = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
    [ Flux ] = psik*JF_star*psi_m.*weight*nuv1;
    Jac_trunc_east_temp = Jac_trunc_east_temp + Flux;

    if data.r_param > 0
        psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
        dspeeddql = 0;%compute_F_dwavespeed_dql(wstar,weast,quadpoint,data.appdata); %???????????????
        F_dql = (JF_star+speedeast.*I_mat - (weast-wstar)*dspeeddql)*0.5;
        %F_dql = compute_F_dflux_dql(wstar,weast,quadpoint,data.appdata);
        dFlux = psik*F_dql*psi_m.*weight*nuv1;
        Jac_east_cell_temp = Jac_east_cell_temp + dFlux;
        psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
        dspeeddqr = 0;%compute_F_dwavespeed_dqr(wstar,weast,quadpoint,data.appdata);        
        F_dqr = (JF_east-speedeast.*I_mat - (weast-wstar)*dspeeddqr)*0.5;  
        %F_dqr = compute_F_dflux_dqr(wstar,weast,quadpoint,data.appdata);
        %????
        dFlux = psik*F_dqr*psi_m.*weight*nuv1;
        Jac_east_other_temp = Jac_east_other_temp + dFlux;
    end
    %west boundary terms
    
    wstar = vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;
    psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    %[speedwest1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    speedwest1 = abs(nuf);
    
    wstar = vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;
    
    %[F_star,JF_star,speedF_star,dspeeddq_star] = problem_compute_all_things(wstar,quadpoint,data.appdata);
    F_star = nuf*wstar; 
    JF_star = nuf; 
    speedF_star = abs(nuf); 
    dspeeddq_star = 0;
    
    
    speedF_west = speedF_star;    
    %[Flux] = problem_F(wstar,quadpoint,data.appdata);
    truncwest_temp = truncwest_temp - nuv1*(psik*F_star)*weight;
    if data.r_param > 0
        %[F_west,JF_west,speedF_west,dspeeddq_west] = problem_compute_all_things(wwest,quadpoint,data.appdata);
        F_west = nuf*wwest; 
        JF_west = nuf; 
        speedF_west = abs(nuf); 
        dspeeddq_west = 0;        
        %[speedF_avg] = problem_F_Jacobianspectralradius((wwest+wstar)*0.5,quadpoint,data.appdata);
        speedF_avg = abs(nuf);
        speedwest = max([speedF_star speedF_west speedF_avg]);
         %[ Flux,speedwest2 ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
        Flux = (F_star+F_west-speedwest.*(wstar-wwest))*0.5;
        fluxwest_temp = fluxwest_temp - nuv1*(psik*Flux)*weight;
    end       
    
    psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
    %F_jac = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
    dFlux = -psik*JF_star*psi_m.*weight*nuv1;
    Jac_trunc_west_temp = Jac_trunc_west_temp + dFlux;

    if data.r_param > 0
        psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
        dspeeddqr = 0;%compute_F_dwavespeed_dqr(wwest,wstar,quadpoint,data.appdata);
        F_dqr = (JF_star-speedwest.*I_mat - (wstar-wwest)*dspeeddqr)*0.5;  
        %F_dqr = compute_F_dflux_dqr(wwest,wstar,quadpoint,data.appdata);
        dFlux = -psik*F_dqr*psi_m.*weight*nuv1;
        Jac_west_cell_temp = Jac_west_cell_temp + dFlux;
        psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
        F_dql = 0;%compute_F_dwavespeed_dql(wwest,wstar,quadpoint,data.appdata);
        dspeeddql = F_dql;
        F_dql = (JF_west-speedwest.*I_mat - (wstar-wwest)*dspeeddql)*0.5;  
        dFlux = -psik*F_dql*psi_m.*weight*nuv1;
        Jac_west_other_temp = Jac_west_other_temp + dFlux;        
    end
    maxspeedF = max([speedF_star speedF_west speedwest1 speedwest2 maxspeedF]);
end

for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        v1quad = data.Dp1list(k,2); 
        weight = data.Dp1quadwgts(k);
        %Derivative terms
        %psim = vectpsi(:,:,tquad,v1quad);                
        wstar = vectpsi(:,:,tquad,v1quad)*qstar;
        
        %[F_star,JF_star,speedF_star,dspeeddq_star] = problem_compute_all_things(wstar,quadpoint,data.appdata);     
        %F_star = data.appdata.nuf*wstar; 
        %JF_star = data.appdata.nuf; 
        %speedF_star = abs(data.appdata.nuf); 
        %dspeeddq_star = 0;
                
        F_star = nuf*wstar; 
        JF_star = nuf; 
        speedF_star = abs(nuf); 
        dspeeddq_star = 0;
        
        
        %psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad);
        %fterm = (nuv1*weight)*psikdxii*F_star;%problem_F(wstar,quadpoint,data.appdata);
        fterm = (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*F_star;%problem_F(wstar,quadpoint,data.appdata);
        
        residual_cell = residual_cell  - fterm;
                
        %F_jac = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
        fterm = weight*nuv1*vectpsi_dxii_Trans(:,:,tquad,v1quad)*JF_star;%F_jac;
        Jac_cell = Jac_cell - fterm*vectpsi(:,:,tquad,v1quad);
        %[speedjac] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        maxspeedF = max(speedF_star,maxspeedF);

end
end