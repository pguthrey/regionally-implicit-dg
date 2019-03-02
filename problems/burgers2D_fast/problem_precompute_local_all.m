function [east_trunc,west_trunc, ...
        Jac_east_trunc,Jac_west_trunc, ...
        nort_trunc,sout_trunc, ...
        Jac_nort_trunc,Jac_sout_trunc, ...
        maxspeedF,maxspeedG] ...
        = problem_precompute_local_all(DGcell,data,~)
% Forms the Jacobian of equations given by the DF data.method for a fiv1ed cell
% written by Pierson Futhrey
% -------------------------------------------------
% INdata.PUTS    DFsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTdata.PUTS   Jacobian
% Note: othr variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

%tauflux = data.tauflux;
%Ipast = data.I_past;
thetaT = data.thetaT;
Pd = data.Pd;
%Pdp1 = data.Pdp1;
%r_param = data.r_param;
Dlist = data.Dlist;
space_dims = data.space_dims;
%Dp1list = data.Dp1list;
%Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv1 = data.nuv1;
nuv2 = data.nuv2;

Jac_east_trunc = zeros(thetaT,thetaT);   
Jac_west_trunc = zeros(thetaT,thetaT); 

east_trunc = zeros(thetaT,1); 
west_trunc = zeros(thetaT,1); 

nort_trunc = zeros(thetaT,1); 
sout_trunc = zeros(thetaT,1); 

Jac_nort_trunc = zeros(thetaT,thetaT);   
Jac_sout_trunc = zeros(thetaT,thetaT); 

maxspeedF = 0;
maxspeedG = 0;

%vectpsi = data.vectpsi;%(:,:,tquad,v1quad,v2quad)
%vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;
%vectpsi_deta_Trans = data.vectpsi_deta_Trans;

vectpsi_west = data.vectpsi_west;%(:,:,tquad,v1quad,v2quad)
vectpsi_east = data.vectpsi_east;%(:,:,tquad,v1quad,v2quad)
%vectpsi_east_Trans = data.vectpsi_east_Trans;
%vectpsi_west_Trans = data.vectpsi_west_Trans;

vectpsi_sout = data.vectpsi_sout;%(:,:,tquad,v1quad,v2quad)
vectpsi_nort = data.vectpsi_nort;%(:,:,tquad,v1quad,v2quad)
%vectpsi_nort_Trans = data.vectpsi_nort_Trans;
%vectpsi_sout_Trans = data.vectpsi_sout_Trans;

q1 = DGcell;

psikpsim_ee_whole = data.psikpsim_ee;
psikpsim_ww_whole = data.psikpsim_ww;

psikpsim_nn_whole = data.psikpsim_nn;
psikpsim_ss_whole = data.psikpsim_ss;


for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    if space_dims >= 2
        v1quad = Dlist(k,2);
    end    
    weight = Dquadwgts(k);    
    wgtnuv1 = weight*nuv1;
    wgtnuv2 = weight*nuv2;
    
    psie = vectpsi_east(:,:,tquad,v1quad);
    psiw = vectpsi_west(:,:,tquad,v1quad);
    psin = vectpsi_nort(:,:,tquad,v1quad);
    psis = vectpsi_sout(:,:,tquad,v1quad);
    
    psieT = psie';
    psiwT = psiw';
    psinT = psin';
    psisT = psis';
    
    psikpsim_ee = psikpsim_ee_whole(:,:,tquad,v1quad);
    psikpsim_ww = psikpsim_ww_whole(:,:,tquad,v1quad);
            
    psikpsim_nn = psikpsim_nn_whole(:,:,tquad,v1quad);
    psikpsim_ss = psikpsim_ss_whole(:,:,tquad,v1quad);   
    
    % EAST AND WEST TERMS
    w1_p = psie*q1;
    w1_m = psiw*q1;
  
    %NORT SOUT terms --------------------------------------------\
    %=============================================================
    w1_n = psin*q1;
    w1_s = psis*q1;
    F1_p = w1_p^2/2;
    F1_m = w1_m^2/2;
    F1_n = w1_n^2/2;
    F1_s = w1_s^2/2;         
    DF1_pDq = w1_p;
    DF1_mDq = w1_m;
    DF1_nDq = w1_n;    
    DF1_sDq = w1_s;
    
    east_trunc = east_trunc + psieT*(F1_p*wgtnuv1);
 
    Jac_east_trunc = Jac_east_trunc + psikpsim_ee*(DF1_pDq*wgtnuv1);
    
    west_trunc = west_trunc - psiwT*(F1_m*wgtnuv1);
    
    Jac_west_trunc = Jac_west_trunc - psikpsim_ww*(DF1_mDq*wgtnuv1); 
           
    nort_trunc = nort_trunc + psinT*(F1_n*wgtnuv2);
    
    Jac_nort_trunc = Jac_nort_trunc + psikpsim_nn*(DF1_nDq*wgtnuv2);
    
    sout_trunc = sout_trunc - psisT*(F1_s*wgtnuv2);
    
    Jac_sout_trunc = Jac_sout_trunc - psikpsim_ss*(DF1_sDq*wgtnuv2);    
end




end



