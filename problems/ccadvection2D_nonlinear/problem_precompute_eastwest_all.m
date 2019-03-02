function [trunceast,truncwest, ...
          fluxeast, fluxwest, ...
          residual_cell, ...
          Jac_east_cell,Jac_west_cell, ...
        Jac_east_other,Jac_west_other, ...
        Jac_trunc_east,Jac_trunc_west, ...
        Jac_cell,maxspeedF] ...
        = problem_precompute_eastwest_all(qstar,qeast,qwest,qpast,data,~)
% Forms the Jacobian of equations given by the DF data.method for a fiv1ed cell
% written by data.Pierson Futhrey
% -------------------------------------------------
% INdata.PUTS    DFsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

residual_cell = data.tauflux*qstar-data.I_past*qpast;
Jac_cell = data.tauflux;

trunceast = zeros(data.thetaT,1);
truncwest = zeros(data.thetaT,1);
fluxeast = zeros(data.thetaT,1);
fluxwest = zeros(data.thetaT,1);

Jac_east_cell = zeros(data.thetaT,data.thetaT);
Jac_west_cell = zeros(data.thetaT,data.thetaT);    
Jac_east_other = zeros(data.thetaT,data.thetaT);
Jac_west_other = zeros(data.thetaT,data.thetaT);  
Jac_trunc_east = zeros(data.thetaT,data.thetaT);   
Jac_trunc_west = zeros(data.thetaT,data.thetaT); 

maxspeedF = 0;

%%{
nuv1 = data.nuv1;

vectpsi = data.vectpsi;
vectpsi_west = data.vectpsi_west;
vectpsi_east = data.vectpsi_east;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west_Trans = data.vectpsi_west_Trans;

Pd = data.Pd;
Pdp1 = data.Pdp1;
r_param = data.r_param;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

v1quad = 1;
v2quad = 1;

for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    if space_dims >= 2
        v1quad = Dlist(k,2);
        if space_dims >= 3
            v2quad = Dlist(k,3);
        end
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
   
    DFwest_pDq = data.appdata.nuf;
    DFeast_mDq = data.appdata.nuf;

    DFstar_pDq = data.appdata.nuf;    
    DFstar_mDq = data.appdata.nuf;

    speedeast = max(abs(data.appdata.nuf),abs(data.appdata.nuf));
    speedwest = max(abs(data.appdata.nuf),abs(data.appdata.nuf));
    
    maxspeedF_quad = max(speedeast,speedwest);
    maxspeedF = max(maxspeedF,maxspeedF_quad);
        
    %east Boundary terms
    psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    
    %trunceast = trunceast + (psik*Fstar_p)*(nuv1*weight);
    if r_param > 0
        Flux = (Fstar_p + Feast_m - speedeast*(weast_m-wstar_p))/2;%  compute_F_numerical_flux(wstar,weast,quadpoint,appdata);
        fluxeast = fluxeast + (psik*Flux)*(nuv1*weight);  
    end
    
    %East Boundary terms
    psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
    Flux = psik*DFstar_pDq*psi_m.*(weight*nuv1);
    Jac_trunc_east = Jac_trunc_east + Flux;
    
    if r_param > 0
        psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
        F_dql = (DFstar_pDq + speedeast)/2;% hmm am I sure ? %compute_F_dflux_dql(wstar,weast,quadpoint,appdata);
        dFlux = psik*F_dql*psi_m.*weight*nuv1;
        Jac_east_cell = Jac_east_cell + dFlux;
        psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
        F_dqr = (DFeast_mDq - speedeast)/2;% compute_F_dflux_dqr(wstar,weast,quadpoint,appdata);
        dFlux = psik*F_dqr*psi_m.*(weight*nuv1);
        Jac_east_other = Jac_east_other + dFlux;
    end    
    
    %west boundary terms
    psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    
    %truncwest = truncwest - (psik*Fstar_m)*(nuv1*weight);
    if r_param > 0
        Flux =  (Fstar_m + Fwest_p - speedwest*(wstar_m-wwest_p))/2;%compute_F_numerical_flux(wwest,wstar,quadpoint,appdata);
        fluxwest = fluxwest - (psik*Flux)*(nuv1*weight);
    end
    psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
    dFlux = psik*DFstar_mDq*psi_m.*(weight*nuv1);
    Jac_trunc_west = Jac_trunc_west - dFlux;
    
    if r_param > 0
        psi_m = vectpsi_west(:,:,tquad,v1quad,v2quad);
        F_dqr = (DFstar_mDq - speedwest)/2;% compute_F_dflux_dqr(wwest,wstar,quadpoint,appdata);
        dFlux = -psik*F_dqr*psi_m.*(weight*nuv1);
        Jac_west_cell = Jac_west_cell + dFlux;
        psi_m = vectpsi_east(:,:,tquad,v1quad,v2quad);
        F_dql = (DFwest_pDq + speedwest)/2;%compute_F_dflux_dql(wwest,wstar,quadpoint,appdata);
        dFlux = -psik*F_dql*psi_m.*weight*nuv1;
        Jac_west_other = Jac_west_other + dFlux;        
    end
    
end

%%{
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
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        Fstar = data.appdata.nuf*wstar;
        gterm = (weight*nuv1)*psikdxii*Fstar;
        residual_cell = residual_cell - gterm; 
        maxspeedF = max(maxspeedF,abs(data.appdata.nuf));
        
        DFstar = data.appdata.nuf;
        psim = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        gterm = (weight*nuv1)*psikdxii*DFstar*psim;
        Jac_cell = Jac_cell - gterm;
end
%}
%{
switch data.M
    case 1
        [trunceast,truncwest,residual_cell, ...
           Jac_trunc_east,Jac_trunc_west,Jac_cell] = problem_exact_integrals_eastwest_M1(qstar,qpast,nuv1,data.appdata.nuf);  
    case 4
        [trunceast,truncwest,residual_cell, ...
           Jac_trunc_east,Jac_trunc_west,Jac_cell] = problem_exact_integrals_eastwest_M4(qstar,qpast,nuv1,data.appdata.nuf);  
    otherwise 
        error('whoops')
end
%}
%}
%disp(['predict eastwest speed = ' num2str(maxspeedF)])



