function [truncnort,truncsout, ...
          fluxnort, fluxsout, ...
          residual_cell, ...
          Jac_nort_cell,Jac_sout_cell, ...
        Jac_nort_other,Jac_sout_other, ...
        Jac_trunc_nort,Jac_trunc_sout, ...
        Jac_cell,maxspeedG] ...
        = problem_precompute_nortsout_all(qstar,qnort,qsout,residual_cell,Jac_cell,data,~)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DFsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

truncnort = zeros(data.thetaT,1); 
truncsout = zeros(data.thetaT,1); 
fluxnort = zeros(data.thetaT,1);
fluxsout = zeros(data.thetaT,1);    

Jac_nort_cell = zeros(data.thetaT,data.thetaT);
Jac_sout_cell = zeros(data.thetaT,data.thetaT);    
Jac_nort_other = zeros(data.thetaT,data.thetaT);
Jac_sout_other = zeros(data.thetaT,data.thetaT);  
Jac_trunc_nort = zeros(data.thetaT,data.thetaT);   
Jac_trunc_sout = zeros(data.thetaT,data.thetaT); 

maxspeedG = 0;

nuv2 = data.nuv2;

vectpsi = data.vectpsi;%(:,:,tquad,v1quad,v2quad)
vectpsi_sout = data.vectpsi_sout;%(:,:,tquad,v1quad,v2quad)
vectpsi_nort = data.vectpsi_nort;%(:,:,tquad,v1quad,v2quad)
vectpsi_deta_Trans = data.vectpsi_deta_Trans;
vectpsi_nort_Trans = data.vectpsi_nort_Trans;
vectpsi_sout_Trans = data.vectpsi_sout_Trans;

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

    Fsout_p = data.appdata.nug*wsout_p;
    Fnort_m = data.appdata.nug*wnort_m;

    Fstar_p = data.appdata.nug*wstar_p;    
    Fstar_m = data.appdata.nug*wstar_m;

    DFsout_pDq = data.appdata.nug;
    DFnort_mDq = data.appdata.nug;

    DFstar_pDq = data.appdata.nug;    
    DFstar_mDq = data.appdata.nug;

    speednort = max(abs(data.appdata.nug),abs(data.appdata.nug));
    speedsout = max(abs(data.appdata.nug),abs(data.appdata.nug));
        
    %nort Boundary terms
    psik = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    truncnort = truncnort + (psik*Fstar_p)*(nuv2*weight);
    
    psi_m = vectpsi_nort(:,:,tquad,v1quad,v2quad);
    Flux = psik*DFstar_pDq*psi_m.*(weight*nuv2);
    Jac_trunc_nort = Jac_trunc_nort + Flux;    
    if r_param > 0
        Flux = (Fstar_p + Fnort_m - speednort*(wnort_m-wstar_p))/2;%  compute_F_numerical_flux(wstar,wnort,quadpoint,appdata);
        fluxnort = fluxnort + (psik*Flux)*(nuv2*weight);  
        
        psi_m = vectpsi_nort(:,:,tquad,v1quad,v2quad);
        F_dql = (DFstar_pDq + speednort)/2;% hmm am I sure ? %compute_F_dflux_dql(wstar,wnort,quadpoint,appdata);
        dFlux = psik*F_dql*psi_m.*weight*nuv2;
        Jac_nort_cell = Jac_nort_cell + dFlux;
        psi_m = vectpsi_sout(:,:,tquad,v1quad,v2quad);
        F_dqr = (DFnort_mDq - speednort)/2;% compute_F_dflux_dqr(wstar,wnort,quadpoint,appdata);
        dFlux = psik*F_dqr*psi_m.*(weight*nuv2);
        Jac_nort_other = Jac_nort_other + dFlux;
    end
    
    %sout boundary terms
    psik = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    truncsout = truncsout - (psik*Fstar_m)*(nuv2*weight);
    
    psi_m = vectpsi_sout(:,:,tquad,v1quad,v2quad);
    dFlux = psik*DFstar_mDq*psi_m.*(weight*nuv2);
    Jac_trunc_sout = Jac_trunc_sout - dFlux;
    if r_param > 0
        Flux =  (Fstar_m + Fsout_p - speedsout*(wstar_m - wsout_p))/2;%compute_F_numerical_flux(wsout,wstar,quadpoint,appdata);
        fluxsout = fluxsout - (psik*Flux)*(nuv2*weight);
        
        psi_m = vectpsi_sout(:,:,tquad,v1quad,v2quad);
        F_dqr = (DFstar_mDq - speedsout)/2;% compute_F_dflux_dqr(wsout,wstar,quadpoint,appdata);
        dFlux = psik*F_dqr*psi_m.*(weight*nuv2);
        Jac_sout_cell = Jac_sout_cell - dFlux;
        psi_m = vectpsi_nort(:,:,tquad,v1quad,v2quad);
        F_dql = (DFsout_pDq + speedsout)/2;%compute_F_dflux_dql(wsout,wstar,quadpoint,appdata);
        dFlux = psik*F_dql*psi_m.*weight*nuv2;
        Jac_sout_other = Jac_sout_other - dFlux;        
    end
    
    maxspeedG_quad = max(speednort,speedsout);
    maxspeedG = max(maxspeedG,maxspeedG_quad);

    maxspeedG_quad = max(speednort,speedsout);
    maxspeedG = max(maxspeedG,maxspeedG_quad);
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
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        Fstar = data.appdata.nug*wstar;
        psikdeta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);

        psim = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        gterm = (weight*nuv2)*psikdeta*Fstar;
        Jac_cell = Jac_cell - gterm*psim;
        maxspeedG = max(maxspeedG,abs(data.appdata.nug));
        gterm = (nuv2*weight)*psikdeta*Fstar;
        residual_cell = residual_cell - gterm; 
        maxspeedG = max(maxspeedG,abs(data.appdata.nug));
end
        %}

%{
switch data.M
    case 1
        [truncnort,truncsout,residual_cell_nortsout, ...
                   Jac_trunc_nort,Jac_trunc_sout,Jac_cell_nortsout] = problem_exact_integrals_nortsout_M1(qstar,nuv2,data.appdata.nug);
    case 4
        [truncnort,truncsout,residual_cell_nortsout, ...
                   Jac_trunc_nort,Jac_trunc_sout,Jac_cell_nortsout] = problem_exact_integrals_nortsout_M4(qstar,nuv2,data.appdata.nug);
    otherwise
        error('whoops')
end
residual_cell = residual_cell + residual_cell_nortsout;
Jac_cell = Jac_cell + Jac_cell_nortsout;
%}

%disp(['predictor nortsout speed = ' num2str(maxspeedG)])
