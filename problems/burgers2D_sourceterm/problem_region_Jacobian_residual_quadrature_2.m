function [Jacobian,residual] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,data)
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


rx = 2;
ry = 2;
q1 = DGregion(:,rx-1,ry-1);
q2 = DGregion(:,rx+0,ry-1);
q3 = DGregion(:,rx+1,ry-1);

q4 = DGregion(:,rx-1,ry+0);
q5 = DGregion(:,rx+0,ry+0);
q6 = DGregion(:,rx+1,ry+0);

q7 = DGregion(:,rx-1,ry+1);
q8 = DGregion(:,rx+0,ry+1);
q9 = DGregion(:,rx+1,ry+1);

q1_past = DGregion_past(:,rx-1,ry-1);
q2_past = DGregion_past(:,rx+0,ry-1);
q3_past = DGregion_past(:,rx+1,ry-1);

q4_past = DGregion_past(:,rx-1,ry+0);
q5_past = DGregion_past(:,rx+0,ry+0);
q6_past = DGregion_past(:,rx+1,ry+0);

q7_past = DGregion_past(:,rx-1,ry+1);
q8_past = DGregion_past(:,rx+0,ry+1);
q9_past = DGregion_past(:,rx+1,ry+1);


if ~data.quadfree_jac_init
    
  %%{
    [Jac_cell_11_test,residual_cell_11_test] = problem_jac_res_init_quadfree_M4_Q(q1,q1_past,nuv1,nuv2);
    [Jac_cell_21_test,residual_cell_21_test] = problem_jac_res_init_quadfree_M4_Q(q2,q2_past,nuv1,nuv2);
    [Jac_cell_31_test,residual_cell_31_test] = problem_jac_res_init_quadfree_M4_Q(q3,q3_past,nuv1,nuv2);
    [Jac_cell_12_test,residual_cell_12_test] = problem_jac_res_init_quadfree_M4_Q(q4,q4_past,nuv1,nuv2);
    [Jac_cell_22_test,residual_cell_22_test] = problem_jac_res_init_quadfree_M4_Q(q5,q5_past,nuv1,nuv2);
    [Jac_cell_32_test,residual_cell_32_test] = problem_jac_res_init_quadfree_M4_Q(q6,q6_past,nuv1,nuv2);
    [Jac_cell_13_test,residual_cell_13_test] = problem_jac_res_init_quadfree_M4_Q(q7,q7_past,nuv1,nuv2);
    [Jac_cell_23_test,residual_cell_23_test] = problem_jac_res_init_quadfree_M4_Q(q8,q8_past,nuv1,nuv2);
    [Jac_cell_33_test,residual_cell_33_test] = problem_jac_res_init_quadfree_M4_Q(q9,q9_past,nuv1,nuv2);
    
    [Jac_cell_11_test] = problem_jac_init_quadfree_M4_Q(q1,nuv1,nuv2);
    [Jac_cell_21_test] = problem_jac_init_quadfree_M4_Q(q2,nuv1,nuv2);
    [Jac_cell_31_test] = problem_jac_init_quadfree_M4_Q(q3,nuv1,nuv2);
    [Jac_cell_12_test] = problem_jac_init_quadfree_M4_Q(q4,nuv1,nuv2);
    [Jac_cell_22_test] = problem_jac_init_quadfree_M4_Q(q5,nuv1,nuv2);
    [Jac_cell_32_test] = problem_jac_init_quadfree_M4_Q(q6,nuv1,nuv2);
    [Jac_cell_13_test] = problem_jac_init_quadfree_M4_Q(q7,nuv1,nuv2);
    [Jac_cell_23_test] = problem_jac_init_quadfree_M4_Q(q8,nuv1,nuv2);
    [Jac_cell_33_test] = problem_jac_init_quadfree_M4_Q(q9,nuv1,nuv2);

    
    %reserror = norm(residual_cell_11_test - residual_cell_11_test);
    %Jacerror = norm(Jac_cell_11_test - Jac_cell_11);

    [residual_cell_11] = problem_res_init(q1,q1_past,data);
    [residual_cell_21] = problem_res_init(q2,q2_past,data);
    [residual_cell_31] = problem_res_init(q3,q3_past,data);
    [residual_cell_12] = problem_res_init(q4,q4_past,data);
    [residual_cell_22] = problem_res_init(q5,q5_past,data);
    [residual_cell_32] = problem_res_init(q6,q6_past,data);
    [residual_cell_13] = problem_res_init(q7,q7_past,data);
    [residual_cell_23] = problem_res_init(q8,q8_past,data);
    [residual_cell_33] = problem_res_init(q9,q9_past,data);
    %}
    [Jac_cell_11,residual_cell_11] = problem_jac_res_init(q1,q1_past,data);
    [Jac_cell_21,residual_cell_21] = problem_jac_res_init(q2,q2_past,data);
    [Jac_cell_31,residual_cell_31] = problem_jac_res_init(q3,q3_past,data);
    [Jac_cell_12,residual_cell_12] = problem_jac_res_init(q4,q4_past,data);
    [Jac_cell_22,residual_cell_22] = problem_jac_res_init(q5,q5_past,data);
    [Jac_cell_32,residual_cell_32] = problem_jac_res_init(q6,q6_past,data);
    [Jac_cell_13,residual_cell_13] = problem_jac_res_init(q7,q7_past,data);
    [Jac_cell_23,residual_cell_23] = problem_jac_res_init(q8,q8_past,data);
    [Jac_cell_33,residual_cell_33] = problem_jac_res_init(q9,q9_past,data);
    
    
    Jacerr = norm(Jac_cell_11 - Jac_cell_11_test);
    reserr = norm(residual_cell_11- residual_cell_11_test);    
    if (Jacerr > 1e-10)||(reserr > 1e-10)
        keyboard 
    end
    Jacerr = norm(Jac_cell_21 - Jac_cell_21_test);
    reserr = norm(residual_cell_21- residual_cell_21_test);    
    if (Jacerr > 1e-10)||(reserr > 1e-10)
        keyboard 
    end
    Jacerr = norm(Jac_cell_31 - Jac_cell_31_test);
    reserr = norm(residual_cell_31- residual_cell_31_test);    
    if (Jacerr > 1e-10)||(reserr > 1e-10)
        keyboard 
    end
    
    Jacerr = norm(Jac_cell_22 - Jac_cell_22_test);
    reserr = norm(residual_cell_22 - residual_cell_22_test);    
    if (Jacerr > 1e-10)||(reserr > 1e-10)
        keyboard 
    end
    
    
else
    switch data.M
        case 4
            [Jac_cell_11,residual_cell_11] = problem_jac_res_init_quadfree_M4_Q(q1,q1_past,nuv1,nuv2);
            [Jac_cell_21,residual_cell_21] = problem_jac_res_init_quadfree_M4_Q(q2,q2_past,nuv1,nuv2);
            [Jac_cell_31,residual_cell_31] = problem_jac_res_init_quadfree_M4_Q(q3,q3_past,nuv1,nuv2);
            [Jac_cell_12,residual_cell_12] = problem_jac_res_init_quadfree_M4_Q(q4,q4_past,nuv1,nuv2);
            [Jac_cell_22,residual_cell_22] = problem_jac_res_init_quadfree_M4_Q(q5,q5_past,nuv1,nuv2);
            [Jac_cell_32,residual_cell_32] = problem_jac_res_init_quadfree_M4_Q(q6,q6_past,nuv1,nuv2);
            [Jac_cell_13,residual_cell_13] = problem_jac_res_init_quadfree_M4_Q(q7,q7_past,nuv1,nuv2);
            [Jac_cell_23,residual_cell_23] = problem_jac_res_init_quadfree_M4_Q(q8,q8_past,nuv1,nuv2);
            [Jac_cell_33,residual_cell_33] = problem_jac_res_init_quadfree_M4_Q(q9,q9_past,nuv1,nuv2);
        case 6
            [Jac_cell_11,residual_cell_11] = problem_jac_res_init_quadfree_M6_Q(q1,q1_past,nuv1,nuv2);
            [Jac_cell_21,residual_cell_21] = problem_jac_res_init_quadfree_M6_Q(q2,q2_past,nuv1,nuv2);
            [Jac_cell_31,residual_cell_31] = problem_jac_res_init_quadfree_M6_Q(q3,q3_past,nuv1,nuv2);
            [Jac_cell_12,residual_cell_12] = problem_jac_res_init_quadfree_M6_Q(q4,q4_past,nuv1,nuv2);
            [Jac_cell_22,residual_cell_22] = problem_jac_res_init_quadfree_M6_Q(q5,q5_past,nuv1,nuv2);
            [Jac_cell_32,residual_cell_32] = problem_jac_res_init_quadfree_M6_Q(q6,q6_past,nuv1,nuv2);
            [Jac_cell_13,residual_cell_13] = problem_jac_res_init_quadfree_M6_Q(q7,q7_past,nuv1,nuv2);
            [Jac_cell_23,residual_cell_23] = problem_jac_res_init_quadfree_M6_Q(q8,q8_past,nuv1,nuv2);
            [Jac_cell_33,residual_cell_33] = problem_jac_res_init_quadfree_M6_Q(q9,q9_past,nuv1,nuv2);
    end
end

psikpsim_ee_whole = data.psikpsim_ee;
psikpsim_we_whole = data.psikpsim_we;
psikpsim_ew_whole = data.psikpsim_ew;
psikpsim_ww_whole = data.psikpsim_ww;

psikpsim_nn_whole = data.psikpsim_nn;
psikpsim_sn_whole = data.psikpsim_sn;
psikpsim_ns_whole = data.psikpsim_ns;
psikpsim_ss_whole = data.psikpsim_ss;

%RESIDUAL
east_flux_23 = 0;
east_flux_22 = 0;
east_flux_21 = 0;
east_flux_13 = 0;
east_flux_11 = 0;
east_flux_12 = 0;

east_trunc_33 = 0;
east_trunc_32 = 0;
east_trunc_31 = 0;

west_flux_33 = 0;
west_flux_32 = 0;
west_flux_31 = 0;
west_flux_23 = 0;
west_flux_22 = 0;
west_flux_21 = 0;


west_trunc_13 = 0;
west_trunc_11 = 0;
west_trunc_12 = 0;

nort_flux_32 = 0;
nort_flux_31 = 0;
nort_flux_22 = 0;
nort_flux_21 = 0;
nort_flux_11 = 0;
nort_flux_12 = 0;


nort_trunc_33 = 0;
nort_trunc_23 = 0;
nort_trunc_13 = 0;

sout_flux_33 = 0;
sout_flux_32 = 0;
sout_flux_23 = 0;
sout_flux_22 = 0;
sout_flux_13 = 0;
sout_flux_12 = 0;


sout_trunc_31 = 0;
sout_trunc_21 = 0;
sout_trunc_11 = 0;

% JACOBIAN 
Jac_east_flux_23 = 0;
Jac_east_flux_22 = 0;
Jac_east_flux_21 = 0;
Jac_east_flux_13 = 0;
Jac_east_flux_11 = 0;
Jac_east_flux_12 = 0;

Jac_east_othr_23 = 0;
Jac_east_othr_22 = 0;
Jac_east_othr_21 = 0;
Jac_east_othr_13 = 0;
Jac_east_othr_11 = 0;
Jac_east_othr_12 = 0;

Jac_east_trunc_33 = 0;
Jac_east_trunc_32 = 0;
Jac_east_trunc_31 = 0;

Jac_west_flux_33 = 0;
Jac_west_flux_32 = 0;
Jac_west_flux_31 = 0;
Jac_west_flux_23 = 0;
Jac_west_flux_22 = 0;
Jac_west_flux_21 = 0;

Jac_west_othr_33 = 0;
Jac_west_othr_32 = 0;
Jac_west_othr_31 = 0;
Jac_west_othr_23 = 0;
Jac_west_othr_22 = 0;
Jac_west_othr_21 = 0;

Jac_west_trunc_13 = 0;
Jac_west_trunc_11 = 0;
Jac_west_trunc_12 = 0;

Jac_nort_flux_32 = 0;
Jac_nort_flux_31 = 0;
Jac_nort_flux_22 = 0;
Jac_nort_flux_21 = 0;
Jac_nort_flux_11 = 0;
Jac_nort_flux_12 = 0;

Jac_nort_othr_32 = 0;
Jac_nort_othr_31 = 0;
Jac_nort_othr_22 = 0;
Jac_nort_othr_21 = 0;
Jac_nort_othr_11 = 0;
Jac_nort_othr_12 = 0;

Jac_nort_trunc_33 = 0;
Jac_nort_trunc_23 = 0;
Jac_nort_trunc_13 = 0;

Jac_sout_flux_33 = 0;
Jac_sout_flux_32 = 0;
Jac_sout_flux_23 = 0;
Jac_sout_flux_22 = 0;
Jac_sout_flux_13 = 0;
Jac_sout_flux_12 = 0;

Jac_sout_othr_33 = 0;
Jac_sout_othr_32 = 0;
Jac_sout_othr_23 = 0;
Jac_sout_othr_22 = 0;
Jac_sout_othr_13 = 0;
Jac_sout_othr_12 = 0;

Jac_sout_trunc_31 = 0;
Jac_sout_trunc_21 = 0;
Jac_sout_trunc_11 = 0;

%{
psie = zeros(1,thetaT);
psiw = zeros(1,thetaT);
psin = zeros(1,thetaT);
psis = zeros(1,thetaT);

psieT = psie';
psiwT = psiw';
psinT = psin';
psisT = psis';

psikpsim_ee = zeros(thetaT,thetaT);
psikpsim_we = zeros(thetaT,thetaT);
psikpsim_ew = zeros(thetaT,thetaT);
psikpsim_ww = zeros(thetaT,thetaT);

psikpsim_nn = zeros(thetaT,thetaT);
psikpsim_sn = zeros(thetaT,thetaT);
psikpsim_ns = zeros(thetaT,thetaT);
psikpsim_ss = zeros(thetaT,thetaT);   
%}

%max = @(l,r) 0 ;

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
    psikpsim_we = psikpsim_we_whole(:,:,tquad,v1quad);
    psikpsim_ew = psikpsim_ew_whole(:,:,tquad,v1quad);
    psikpsim_ww = psikpsim_ww_whole(:,:,tquad,v1quad);
            
    psikpsim_nn = psikpsim_nn_whole(:,:,tquad,v1quad);
    psikpsim_sn = psikpsim_sn_whole(:,:,tquad,v1quad);
    psikpsim_ns = psikpsim_ns_whole(:,:,tquad,v1quad);
    psikpsim_ss = psikpsim_ss_whole(:,:,tquad,v1quad);   
    
    % EAST AND WEST TERMS
    w1_p = psie*q1;
    w1_m = psiw*q1;
    
    w2_p = psie*q2;
    w2_m = psiw*q2;
    
    w3_p = psie*q3;
    w3_m = psiw*q3;    
    
    w4_p = psie*q4;
    w4_m = psiw*q4;    

    w5_p = psie*q5;
    w5_m = psiw*q5;    

    w6_p = psie*q6;
    w6_m = psiw*q6;    

    w7_p = psie*q7;
    w7_m = psiw*q7;    

    w8_p = psie*q8;
    w8_m = psiw*q8;    

    w9_p = psie*q9;
    w9_m = psiw*q9;    
    

    %NORT SOUT terms --------------------------------------------\
    %=============================================================
    w1_n = psin*q1;
    w1_s = psis*q1;

    w2_n = psin*q2;
    w2_s = psis*q2;    

    w3_n = psin*q3;
    w3_s = psis*q3;    

    w4_n = psin*q4;
    w4_s = psis*q4;    

    w5_n = psin*q5;
    w5_s = psis*q5;    

    w6_n = psin*q6;
    w6_s = psis*q6;    

    w7_n = psin*q7;
    w7_s = psis*q7;    

    w8_n = psin*q8;
    w8_s = psis*q8;    

    w9_n = psin*q9;
    w9_s = psis*q9;    
    
    F1_p = w1_p^2/2;
    F2_p = w2_p^2/2;
    F3_p = w3_p^2/2;
    F4_p = w4_p^2/2;
    F5_p = w5_p^2/2;
    F6_p = w6_p^2/2;
    F7_p = w7_p^2/2;
    F8_p = w8_p^2/2;
    F9_p = w9_p^2/2;

    F1_m = w1_m^2/2;
    F2_m = w2_m^2/2;
    F3_m = w3_m^2/2;
    F4_m = w4_m^2/2;
    F5_m = w5_m^2/2;
    F6_m = w6_m^2/2;
    F7_m = w7_m^2/2;
    F8_m = w8_m^2/2;
    F9_m = w9_m^2/2;

    F1_n = w1_n^2/2;
    F2_n = w2_n^2/2;
    F3_n = w3_n^2/2;
    F4_n = w4_n^2/2;
    F5_n = w5_n^2/2;
    F6_n = w6_n^2/2;
    F7_n = w7_n^2/2;
    F8_n = w8_n^2/2;
    F9_n = w9_n^2/2;

    F1_s = w1_s^2/2;
    F2_s = w2_s^2/2;
    F3_s = w3_s^2/2;
    F4_s = w4_s^2/2;
    F5_s = w5_s^2/2;
    F6_s = w6_s^2/2;
    F7_s = w7_s^2/2;
    F8_s = w8_s^2/2;
    F9_s = w9_s^2/2;
         
    DF1_pDq = w1_p;
    DF2_pDq = w2_p;
    DF3_pDq = w3_p;
    DF4_pDq = w4_p;
    DF5_pDq = w5_p;
    DF6_pDq = w6_p;
    DF7_pDq = w7_p;
    DF8_pDq = w8_p;
    DF9_pDq = w9_p;
    
    DF1_mDq = w1_m;
    DF2_mDq = w2_m;
    DF3_mDq = w3_m;
    DF4_mDq = w4_m;
    DF5_mDq = w5_m;
    DF6_mDq = w6_m;
    DF7_mDq = w7_m;
    DF8_mDq = w8_m;
    DF9_mDq = w9_m;    
    
    DF1_nDq = w1_n;
    DF2_nDq = w2_n;
    DF3_nDq = w3_n;
    DF4_nDq = w4_n;
    DF5_nDq = w5_n;
    DF6_nDq = w6_n;
    DF7_nDq = w7_n;
    DF8_nDq = w8_n;
    DF9_nDq = w9_n;
    
    DF1_sDq = w1_s;
    DF2_sDq = w2_s;
    DF3_sDq = w3_s;
    DF4_sDq = w4_s;
    DF5_sDq = w5_s;
    DF6_sDq = w6_s;
    DF7_sDq = w7_s;
    DF8_sDq = w8_s;
    DF9_sDq = w9_s;   
    
    %psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    east_trunc_31 = east_trunc_31 + psieT*(F3_p*wgtnuv1);
    east_trunc_32 = east_trunc_32 + psieT*(F6_p*wgtnuv1);
    east_trunc_33 = east_trunc_33 + psieT*(F9_p*wgtnuv1);

    Jac_east_trunc_31 = Jac_east_trunc_31 + psikpsim_ee*(DF3_pDq*wgtnuv1);
    Jac_east_trunc_32 = Jac_east_trunc_32 + psikpsim_ee*(DF6_pDq*wgtnuv1);
    Jac_east_trunc_33 = Jac_east_trunc_33 + psikpsim_ee*(DF9_pDq*wgtnuv1);      
    
    west_trunc_11 = west_trunc_11 - psiwT*(F1_m*wgtnuv1);
    west_trunc_12 = west_trunc_12 - psiwT*(F4_m*wgtnuv1);
    west_trunc_13 = west_trunc_13 - psiwT*(F7_m*wgtnuv1);
    
    Jac_west_trunc_11 = Jac_west_trunc_11 - psikpsim_ww*(DF1_mDq*wgtnuv1);
    Jac_west_trunc_12 = Jac_west_trunc_12 - psikpsim_ww*(DF4_mDq*wgtnuv1);
    Jac_west_trunc_13 = Jac_west_trunc_13 - psikpsim_ww*(DF7_mDq*wgtnuv1);
          
    % East terms
    % Interface 1-2
    speed = max(abs(w1_p),abs(w2_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F1_p + F2_m - speed*(w2_m-w1_p))*wgtnuv1/2;
    F_dql = (DF1_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF2_mDq - speed)*wgtnuv1/2;    
    east_flux_11 = east_flux_11 + psieT*Flux;
    west_flux_21 = west_flux_21 - psiwT*Flux;
    Jac_east_flux_11 = Jac_east_flux_11 + psikpsim_ee*F_dql; 
    Jac_east_othr_11 = Jac_east_othr_11 + psikpsim_ew*F_dqr;    
    Jac_west_flux_21 = Jac_west_flux_21 - psikpsim_ww*F_dqr; 
    Jac_west_othr_21 = Jac_west_othr_21 - psikpsim_we*F_dql; 
        
    % Interface 2-3
    speed = max(abs(w2_p),abs(w3_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F2_p + F3_m - speed*(w3_m-w2_p))*wgtnuv1/2;
    F_dql = (DF2_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF3_mDq - speed)*wgtnuv1/2;
    east_flux_21 = east_flux_21 + psieT*Flux;
    west_flux_31 = west_flux_31 - psiwT*Flux;
    Jac_east_flux_21 = Jac_east_flux_21 + psikpsim_ee*F_dql; 
    Jac_east_othr_21 = Jac_east_othr_21 + psikpsim_ew*F_dqr; 
    Jac_west_flux_31 = Jac_west_flux_31 - psikpsim_ww*F_dqr; 
    Jac_west_othr_31 = Jac_west_othr_31 - psikpsim_we*F_dql; 
    
    % Interface 4-5
    speed = max(abs(w4_p),abs(w5_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F4_p + F5_m - speed*(w5_m-w4_p))*wgtnuv1/2;
    F_dql = (DF4_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF5_mDq - speed)*wgtnuv1/2;
    east_flux_12 = east_flux_12 + psieT*Flux;
    west_flux_22 = west_flux_22 - psiwT*Flux;
    Jac_east_flux_12 = Jac_east_flux_12 + psikpsim_ee*F_dql; 
    Jac_east_othr_12 = Jac_east_othr_12 + psikpsim_ew*F_dqr; 
    Jac_west_flux_22 = Jac_west_flux_22 - psikpsim_ww*F_dqr; 
    Jac_west_othr_22 = Jac_west_othr_22 - psikpsim_we*F_dql; 
    
    
    % Interface 5-6
    speed = max(abs(w5_p),abs(w6_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F5_p + F6_m - speed*(w6_m-w5_p))*wgtnuv1/2;
    east_flux_22 = east_flux_22 + psieT*Flux;
    west_flux_32 = west_flux_32 - psiwT*Flux;
    F_dql = (DF5_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF6_mDq - speed)*wgtnuv1/2;
    Jac_east_flux_22 = Jac_east_flux_22 + psikpsim_ee*F_dql; 
    Jac_east_othr_22 = Jac_east_othr_22 + psikpsim_ew*F_dqr; 
    Jac_west_flux_32 = Jac_west_flux_32 - psikpsim_ww*F_dqr; 
    Jac_west_othr_32 = Jac_west_othr_32 - psikpsim_we*F_dql; 

    
    
    % Interface 7-8
    speed = max(abs(w7_p),abs(w8_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F7_p + F8_m - speed*(w8_m-w7_p))*wgtnuv1/2;
    F_dql = (DF7_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF8_mDq - speed)*wgtnuv1/2;
    east_flux_13 = east_flux_13 + psieT*Flux;
    west_flux_23 = west_flux_23 - psiwT*Flux;
    Jac_east_flux_13 = Jac_east_flux_13 + psikpsim_ee*F_dql; 
    Jac_east_othr_13 = Jac_east_othr_13 + psikpsim_ew*F_dqr; 
    Jac_west_flux_23 = Jac_west_flux_23 - psikpsim_ww*F_dqr; 
    Jac_west_othr_23 = Jac_west_othr_23 - psikpsim_we*F_dql;    
    
    
    % Interface 8-9
    speed = max(abs(w8_p),abs(w9_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F8_p + F9_m - speed*(w9_m-w8_p))*wgtnuv1/2;
    F_dql = (DF8_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF9_mDq - speed)*wgtnuv1/2;
    east_flux_23 = east_flux_23 + psieT*Flux;
    west_flux_33 = west_flux_33 - psiwT*Flux;
    Jac_east_flux_23 = Jac_east_flux_23 + psikpsim_ee*F_dql; 
    Jac_east_othr_23 = Jac_east_othr_23 + psikpsim_ew*F_dqr;    
    Jac_west_flux_33 = Jac_west_flux_33 - psikpsim_ww*F_dqr; 
    Jac_west_othr_33 = Jac_west_othr_33 - psikpsim_we*F_dql; 
        
    nort_trunc_13 = nort_trunc_13 + psinT*(F7_n*wgtnuv2);
    nort_trunc_23 = nort_trunc_23 + psinT*(F8_n*wgtnuv2);
    nort_trunc_33 = nort_trunc_33 + psinT*(F9_n*wgtnuv2);
    
    Jac_nort_trunc_13 = Jac_nort_trunc_13 + psikpsim_nn*(DF7_nDq*wgtnuv2);
    Jac_nort_trunc_23 = Jac_nort_trunc_23 + psikpsim_nn*(DF8_nDq*wgtnuv2);
    Jac_nort_trunc_33 = Jac_nort_trunc_33 + psikpsim_nn*(DF9_nDq*wgtnuv2);      
    
    sout_trunc_11 = sout_trunc_11 - psisT*(F1_s*wgtnuv2);
    sout_trunc_21 = sout_trunc_21 - psisT*(F2_s*wgtnuv2);
    sout_trunc_31 = sout_trunc_31 - psisT*(F3_s*wgtnuv2);
    
    Jac_sout_trunc_11 = Jac_sout_trunc_11 - psikpsim_ss*(DF1_sDq*wgtnuv2);
    Jac_sout_trunc_21 = Jac_sout_trunc_21 - psikpsim_ss*(DF2_sDq*wgtnuv2);
    Jac_sout_trunc_31 = Jac_sout_trunc_31 - psikpsim_ss*(DF3_sDq*wgtnuv2);
    
    % nort terms -------------------------------------------------------
    % Interface 1-4
    speed = max(abs(w1_n),abs(w4_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F1_n + F4_s - speed*(w4_s-w1_n))*wgtnuv2/2;
    F_dql = (DF1_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF4_sDq - speed)*wgtnuv2/2;
    nort_flux_11 = nort_flux_11 + psinT*Flux;
    sout_flux_12 = sout_flux_12 - psisT*Flux;
    Jac_nort_flux_11 = Jac_nort_flux_11 + psikpsim_nn*F_dql;    
    Jac_nort_othr_11 = Jac_nort_othr_11 + psikpsim_ns*F_dqr;    
    Jac_sout_flux_12 = Jac_sout_flux_12 - psikpsim_ss*F_dqr;    
    Jac_sout_othr_12 = Jac_sout_othr_12 - psikpsim_sn*F_dql;    
        
    
    % Interface 2-5
    speed = max(abs(w2_n),abs(w5_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F2_n + F5_s - speed*(w5_s-w2_n))*wgtnuv2/2;
    F_dql = (DF2_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF5_sDq - speed)*wgtnuv2/2;
    nort_flux_21 = nort_flux_21 + psinT*Flux;
    sout_flux_22 = sout_flux_22 - psisT*Flux;
    Jac_nort_flux_21 = Jac_nort_flux_21 + psikpsim_nn*F_dql;    
    Jac_nort_othr_21 = Jac_nort_othr_21 + psikpsim_ns*F_dqr;      
    Jac_sout_flux_22 = Jac_sout_flux_22 - psikpsim_ss*F_dqr;    
    Jac_sout_othr_22 = Jac_sout_othr_22 - psikpsim_sn*F_dql;  
        
    
    % Interface 3-6
    speed = max(abs(w3_n),abs(w6_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F3_n + F6_s - speed*(w6_s-w3_n))*wgtnuv2/2;
    F_dql = (DF3_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF6_sDq - speed)*wgtnuv2/2;
    nort_flux_31 = nort_flux_31 + psinT*Flux;
    sout_flux_32 = sout_flux_32 - psisT*Flux;
    Jac_nort_flux_31 = Jac_nort_flux_31 + psikpsim_nn*F_dql;    
    Jac_nort_othr_31 = Jac_nort_othr_31 + psikpsim_ns*F_dqr;     
    Jac_sout_flux_32 = Jac_sout_flux_32 - psikpsim_ss*F_dqr;    
    Jac_sout_othr_32 = Jac_sout_othr_32 - psikpsim_sn*F_dql;    
    
    
    % Interface 4-7
    speed = max(abs(w4_n),abs(w7_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F4_n + F7_s - speed*(w7_s-w4_n))*wgtnuv2/2;
    F_dql = (DF4_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF7_sDq - speed)*wgtnuv2/2;
    nort_flux_12 = nort_flux_12 + psinT*Flux;
    sout_flux_13 = sout_flux_13 - psisT*Flux;
    Jac_nort_flux_12 = Jac_nort_flux_12 + psikpsim_nn*F_dql;    
    Jac_nort_othr_12 = Jac_nort_othr_12 + psikpsim_ns*F_dqr;    
    Jac_sout_flux_13 = Jac_sout_flux_13 - psikpsim_ss*F_dqr;    
    Jac_sout_othr_13 = Jac_sout_othr_13 - psikpsim_sn*F_dql;  
    
    
    % Interface 5-8
    speed = max(abs(w5_n),abs(w8_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F5_n + F8_s - speed*(w8_s-w5_n))*wgtnuv2/2;
    F_dql = (DF5_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF8_sDq - speed)*wgtnuv2/2;
    nort_flux_22 = nort_flux_22 + psinT*Flux;
    sout_flux_23 = sout_flux_23 - psisT*Flux;
    Jac_nort_flux_22 = Jac_nort_flux_22 + psikpsim_nn*F_dql;    
    Jac_nort_othr_22 = Jac_nort_othr_22 + psikpsim_ns*F_dqr;     
    Jac_sout_flux_23 = Jac_sout_flux_23 - psikpsim_ss*F_dqr;    
    Jac_sout_othr_23 = Jac_sout_othr_23 - psikpsim_sn*F_dql;     
    
    % Interface 6-9
    speed = max(abs(w6_n),abs(w9_s));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F6_n + F9_s - speed*(w9_s-w6_n))*wgtnuv2/2;
    nort_flux_32 = nort_flux_32 + psinT*Flux;
    sout_flux_33 = sout_flux_33 - psisT*Flux;
    F_dql = (DF6_nDq + speed)*wgtnuv2/2;
    F_dqr = (DF9_sDq - speed)*wgtnuv2/2;
    Jac_nort_flux_32 = Jac_nort_flux_32 + psikpsim_nn*F_dql;
    Jac_nort_othr_32 = Jac_nort_othr_32 + psikpsim_ns*F_dqr;
    Jac_sout_flux_33 = Jac_sout_flux_33 - psikpsim_ss*F_dqr;
    Jac_sout_othr_33 = Jac_sout_othr_33 - psikpsim_sn*F_dql;
    
end

indices = @(i) (1:thetaT) + (i-1)*thetaT;
                               
residual(indices(1),1) =   residual_cell_11 + ...
                                    east_flux_11 + ... 
                                    west_trunc_11 + ...
                                    nort_flux_11 + ... 
                                    sout_trunc_11 ;
                                
residual(indices(2),1) =   residual_cell_21 + ...
                                    east_flux_21 + ... 
                                    west_flux_21 + ...
                                    nort_flux_21 + ... 
                                    sout_trunc_21 ;
residual(indices(3),1) =   residual_cell_31 + ...
                                    east_trunc_31 + ... 
                                    west_flux_31 + ...
                                    nort_flux_31 + ... 
                                    sout_trunc_31 ;
                                
residual(indices(4),1) =   residual_cell_12 + ...
                                    east_flux_12 + ... 
                                    west_trunc_12 + ...
                                    nort_flux_12 + ... 
                                    sout_flux_12 ;
residual(indices(5),1) =   residual_cell_22 + ...
                                    east_flux_22 + ... 
                                    west_flux_22 + ...
                                    nort_flux_22 + ... 
                                    sout_flux_22 ;
residual(indices(6),1) =   residual_cell_32 + ...
                                    east_trunc_32 + ... 
                                    west_flux_32 + ...
                                    nort_flux_32 + ... 
                                    sout_flux_32 ;

residual(indices(7),1) =   residual_cell_13 + ...
                                    east_flux_13 + ... 
                                    west_trunc_13 + ...
                                    nort_trunc_13 + ... 
                                    sout_flux_13 ;
                                
residual(indices(8),1) =   residual_cell_23 + ...
                                    east_flux_23 + ... 
                                    west_flux_23 + ...
                                    nort_trunc_23 + ... 
                                    sout_flux_23 ;
residual(indices(9),1) =   residual_cell_33 + ...
                                    east_trunc_33 + ... 
                                    west_flux_33 + ...
                                    nort_trunc_33 + ... 
                                    sout_flux_33 ;


Jacobian(indices(1),indices(1)) =   Jac_cell_11 + ...
                                    Jac_east_flux_11 + ... 
                                    Jac_west_trunc_11 + ...
                                    Jac_nort_flux_11 + ... 
                                    Jac_sout_trunc_11 ;
Jacobian(indices(1),indices(2)) = Jac_east_othr_11;
Jacobian(indices(1),indices(4)) = Jac_nort_othr_11;

                                
                                
Jacobian(indices(2),indices(2)) =   Jac_cell_21 + ...
                                    Jac_east_flux_21 + ... 
                                    Jac_west_flux_21 + ...
                                    Jac_nort_flux_21 + ... 
                                    Jac_sout_trunc_21 ;
Jacobian(indices(2),indices(1)) = Jac_west_othr_21;
Jacobian(indices(2),indices(3)) = Jac_east_othr_21;
Jacobian(indices(2),indices(5)) = Jac_nort_othr_21;
                                
                                
                                
Jacobian(indices(3),indices(3)) =   Jac_cell_31 + ...
                                    Jac_east_trunc_31 + ... 
                                    Jac_west_flux_31 + ...
                                    Jac_nort_flux_31 + ... 
                                    Jac_sout_trunc_31 ;
Jacobian(indices(3),indices(2)) = Jac_west_othr_31;
Jacobian(indices(3),indices(6)) = Jac_nort_othr_31;
                                
                                
                                
Jacobian(indices(4),indices(4)) =   Jac_cell_12 + ...
                                    Jac_east_flux_12 + ... 
                                    Jac_west_trunc_12 + ...
                                    Jac_nort_flux_12 + ... 
                                    Jac_sout_flux_12 ;
Jacobian(indices(4),indices(1)) = Jac_sout_othr_12;
Jacobian(indices(4),indices(5)) = Jac_east_othr_12;
Jacobian(indices(4),indices(7)) = Jac_nort_othr_12;
                                
                                
Jacobian(indices(5),indices(5)) =   Jac_cell_22 + ...
                                    Jac_east_flux_22 + ... 
                                    Jac_west_flux_22 + ...
                                    Jac_nort_flux_22 + ... 
                                    Jac_sout_flux_22 ;
Jacobian(indices(5),indices(2)) = Jac_sout_othr_22;
Jacobian(indices(5),indices(4)) = Jac_west_othr_22;
Jacobian(indices(5),indices(6)) = Jac_east_othr_22;
Jacobian(indices(5),indices(8)) = Jac_nort_othr_22;

                                
                                
Jacobian(indices(6),indices(6)) =   Jac_cell_32 + ...
                                    Jac_east_trunc_32 + ... 
                                    Jac_west_flux_32 + ...
                                    Jac_nort_flux_32 + ... 
                                    Jac_sout_flux_32 ;
Jacobian(indices(6),indices(3)) = Jac_sout_othr_32;
Jacobian(indices(6),indices(5)) = Jac_west_othr_32;
Jacobian(indices(6),indices(9)) = Jac_nort_othr_32;
                                

Jacobian(indices(7),indices(7)) =   Jac_cell_13 + ...
                                    Jac_east_flux_13 + ... 
                                    Jac_west_trunc_13 + ...
                                    Jac_nort_trunc_13 + ... 
                                    Jac_sout_flux_13 ;
Jacobian(indices(7),indices(4)) = Jac_sout_othr_13;
Jacobian(indices(7),indices(8)) = Jac_east_othr_13;
                                
                                
Jacobian(indices(8),indices(8)) =   Jac_cell_23 + ...
                                    Jac_east_flux_23 + ... 
                                    Jac_west_flux_23 + ...
                                    Jac_nort_trunc_23 + ... 
                                    Jac_sout_flux_23 ;
Jacobian(indices(8),indices(5)) = Jac_sout_othr_23;
Jacobian(indices(8),indices(7)) = Jac_west_othr_23;
Jacobian(indices(8),indices(9)) = Jac_east_othr_23;
                                
                                
                                
Jacobian(indices(9),indices(9)) =   Jac_cell_33 + ...
                                    Jac_east_trunc_33 + ... 
                                    Jac_west_flux_33 + ...
                                    Jac_nort_trunc_33 + ... 
                                    Jac_sout_flux_33 ;
Jacobian(indices(9),indices(6)) = Jac_sout_othr_33;
Jacobian(indices(9),indices(8)) = Jac_west_othr_33;
                                



end



