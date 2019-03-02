function [Jac_east_cell_temp,Jac_west_cell_temp, ...
        Jac_east_other_temp,Jac_west_other_temp, ...
        Jac_trunc_east_temp,Jac_trunc_west_temp, ...
        Jac_cell,maxspeedF] ...
        = predictor_precompute_eastwest_Jacobian(qstar,qeast,qwest,data,cellcenter)
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
    [speedeast1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);

    psi_m = data.vectpsi_east(:,:,tquad,v1quad,v2quad);
    F_jac = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
    [ Flux ] = psik*F_jac*psi_m.*weight*data.nuv1;
    Jac_trunc_east_temp = Jac_trunc_east_temp + Flux;

    if data.r_param > 0
        psi_m = data.vectpsi_east(:,:,tquad,v1quad,v2quad);
        F_dql = compute_F_dflux_dql(wstar,weast,quadpoint,data.appdata);
        dFlux = psik*F_dql*psi_m.*weight*data.nuv1;
        Jac_east_cell_temp = Jac_east_cell_temp + dFlux;
        psi_m = data.vectpsi_west(:,:,tquad,v1quad,v2quad);
        F_dqr = compute_F_dflux_dqr(wstar,weast,quadpoint,data.appdata);
        dFlux = psik*F_dqr*psi_m.*weight*data.nuv1;
        Jac_east_other_temp = Jac_east_other_temp + dFlux;
    end
    %west boundary terms
    v1 = v1center-data.deltav1/2;        
    v2 = v2center+data.deltav2/2*v1loc;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];

    wstar = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    [speedwest1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    
    psi_m = data.vectpsi_west(:,:,tquad,v1quad,v2quad);
    F_jac = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
    dFlux = -psik*F_jac*psi_m.*weight*data.nuv1;
    Jac_trunc_west_temp = Jac_trunc_west_temp + dFlux;

    if data.r_param > 0
        psi_m = data.vectpsi_west(:,:,tquad,v1quad,v2quad);
        F_dqr = compute_F_dflux_dqr(wwest,wstar,quadpoint,data.appdata);
        dFlux = -psik*F_dqr*psi_m.*weight*data.nuv1;
        Jac_west_cell_temp = Jac_west_cell_temp + dFlux;
        psi_m = data.vectpsi_east(:,:,tquad,v1quad,v2quad);
        F_dql = compute_F_dflux_dql(wwest,wstar,quadpoint,data.appdata);
        dFlux = -psik*F_dql*psi_m.*weight*data.nuv1;
        Jac_west_other_temp = Jac_west_other_temp + dFlux;        
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
        psim = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        F_jac = problem_F_FluxJacobian(wstar,quadpoint,data.appdata);
        fterm = weight*data.nuv1*psikdxii*F_jac;
        Jac_cell = Jac_cell - fterm*psim;
        [speedjac] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        maxspeedF = max([speedjac maxspeedF]);

end
end