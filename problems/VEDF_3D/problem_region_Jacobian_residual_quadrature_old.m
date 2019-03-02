function [Jacobian,residual] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_east,data)
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
%Ipast = data.I_east;
thetaT = data.thetaT;
Pd = data.Pd;
%Pdp1 = data.Pdp1;
%r_earam = data.r_earam;
Dlist = data.Dlist;
space_dims = data.space_dims;
%Dp1list = data.Dp1list;
%Dp1quadwgts = data.Dp1quadwgts;
uadwgts = data.uadwgts;

nuv1 = data.nuv1;
nuv2 = data.nuv2;
nuv3 = data.nuv3;

maxspeedF = 0;
maxspeedG = 0;
maxspeedH = 0;

vectpsi = data.vectpsi;

vectpsi_west = data.vectpsi_west;
vectpsi_east = data.vectpsi_east;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west_Trans = data.vectpsi_west_Trans;

vectpsi_sout = data.vectpsi_sout;
vectpsi_nort = data.vectpsi_nort;
vectpsi_nort_Trans = data.vectpsi_nort_Trans;
vectpsi_sout_Trans = data.vectpsi_sout_Trans;

vectpsi_down = data.vectpsi_down;
vectpsi_uppr = data.vectpsi_uppr;
vectpsi_uppr_Trans = data.vectpsi_uppr_Trans;
vectpsi_down_Trans = data.vectpsi_down_Trans;

rx = 2;
ry = 2;
rz = 2;

q01 = DGregion(:,rx-1,ry-1,rz-1);
q02 = DGregion(:,rx+0,ry-1,rz-1);
q03 = DGregion(:,rx+1,ry-1,rz-1);

q04 = DGregion(:,rx-1,ry+0,rz-1);
q05 = DGregion(:,rx+0,ry+0,rz-1);
q06 = DGregion(:,rx+1,ry+0,rz-1);

q07 = DGregion(:,rx-1,ry+1,rz-1);
q08 = DGregion(:,rx+0,ry+1,rz-1);
q09 = DGregion(:,rx+1,ry+1,rz-1);

q10 = DGregion(:,rx-1,ry-1,rz);
q11 = DGregion(:,rx+0,ry-1,rz);
q12 = DGregion(:,rx+1,ry-1,rz);

q13 = DGregion(:,rx-1,ry+0,rz);
q14 = DGregion(:,rx+0,ry+0,rz);
q15 = DGregion(:,rx+1,ry+0,rz);

q16 = DGregion(:,rx-1,ry+1,rz);
q17 = DGregion(:,rx+0,ry+1,rz);
q18 = DGregion(:,rx+1,ry+1,rz);

q19 = DGregion(:,rx-1,ry-1,rz+1);
q20 = DGregion(:,rx+0,ry-1,rz+1);
q21 = DGregion(:,rx+1,ry-1,rz+1);

q22 = DGregion(:,rx-1,ry+0,rz+1);
q23 = DGregion(:,rx+0,ry+0,rz+1);
q24 = DGregion(:,rx+1,ry+0,rz+1);

q25 = DGregion(:,rx-1,ry+1,rz+1);
q26 = DGregion(:,rx+0,ry+1,rz+1);
q27 = DGregion(:,rx+1,ry+1,rz+1);

qpast01 = DGregion_east(:,rx-1,ry-1,rz-1);
qpast02 = DGregion_east(:,rx+0,ry-1,rz-1);
qpast03 = DGregion_east(:,rx+1,ry-1,rz-1);

qpast04 = DGregion_east(:,rx-1,ry+0,rz-1);
qpast05 = DGregion_east(:,rx+0,ry+0,rz-1);
qpast06 = DGregion_east(:,rx+1,ry+0,rz-1);

qpast07 = DGregion_east(:,rx-1,ry+1,rz-1);
qpast08 = DGregion_east(:,rx+0,ry+1,rz-1);
qpast09 = DGregion_east(:,rx+1,ry+1,rz-1);

qpast10 = DGregion_east(:,rx-1,ry-1,rz);
qpast11 = DGregion_east(:,rx+0,ry-1,rz);
qpast12 = DGregion_east(:,rx+1,ry-1,rz);

qpast13 = DGregion_east(:,rx-1,ry+0,rz);
qpast14 = DGregion_east(:,rx+0,ry+0,rz);
qpast15 = DGregion_east(:,rx+1,ry+0,rz);

qpast16 = DGregion_east(:,rx-1,ry+1,rz);
qpast17 = DGregion_east(:,rx+0,ry+1,rz);
qpast18 = DGregion_east(:,rx+1,ry+1,rz);

qpast19 = DGregion_east(:,rx-1,ry-1,rz+1);
qpast20 = DGregion_east(:,rx+0,ry-1,rz+1);
qpast21 = DGregion_east(:,rx+1,ry-1,rz+1);

qpast22 = DGregion_east(:,rx-1,ry+0,rz+1);
qpast23 = DGregion_east(:,rx+0,ry+0,rz+1);
qpast24 = DGregion_east(:,rx+1,ry+0,rz+1);

qpast25 = DGregion_east(:,rx-1,ry+1,rz+1);
qpast26 = DGregion_east(:,rx+0,ry+1,rz+1);
qpast27 = DGregion_east(:,rx+1,ry+1,rz+1);

[Jac_cell_111,residual_cell_111] = problem_jac_res_init(q01,q01_east,data);
[Jac_cell_211,residual_cell_211] = problem_jac_res_init(q02,q02_east,data);
[Jac_cell_311,residual_cell_311] = problem_jac_res_init(q03,q03_east,data);
[Jac_cell_121,residual_cell_121] = problem_jac_res_init(q04,q04_east,data);
[Jac_cell_221,residual_cell_221] = problem_jac_res_init(q05,q05_east,data);
[Jac_cell_321,residual_cell_321] = problem_jac_res_init(q06,q06_east,data);
[Jac_cell_131,residual_cell_131] = problem_jac_res_init(q07,q07_east,data);
[Jac_cell_231,residual_cell_231] = problem_jac_res_init(q08,q08_east,data);
[Jac_cell_331,residual_cell_331] = problem_jac_res_init(q09,q09_east,data);

[Jac_cell_112,residual_cell_112] = problem_jac_res_init(q10,q10_east,data);
[Jac_cell_212,residual_cell_212] = problem_jac_res_init(q11,q11_east,data);
[Jac_cell_312,residual_cell_312] = problem_jac_res_init(q12,q12_east,data);
[Jac_cell_122,residual_cell_122] = problem_jac_res_init(q13,q13_east,data);
[Jac_cell_222,residual_cell_222] = problem_jac_res_init(q14,q14_east,data);
[Jac_cell_322,residual_cell_322] = problem_jac_res_init(q15,q15_east,data);
[Jac_cell_132,residual_cell_132] = problem_jac_res_init(q16,q16_east,data);
[Jac_cell_232,residual_cell_232] = problem_jac_res_init(q17,q17_east,data);
[Jac_cell_332,residual_cell_332] = problem_jac_res_init(q18,q18_east,data);

[Jac_cell_113,residual_cell_113] = problem_jac_res_init(q19,q19_east,data);
[Jac_cell_213,residual_cell_213] = problem_jac_res_init(q20,q20_east,data);
[Jac_cell_313,residual_cell_313] = problem_jac_res_init(q21,q21_east,data);
[Jac_cell_123,residual_cell_123] = problem_jac_res_init(q22,q22_east,data);
[Jac_cell_223,residual_cell_223] = problem_jac_res_init(q23,q23_east,data);
[Jac_cell_323,residual_cell_323] = problem_jac_res_init(q24,q24_east,data);
[Jac_cell_133,residual_cell_133] = problem_jac_res_init(q25,q25_east,data);
[Jac_cell_233,residual_cell_233] = problem_jac_res_init(q26,q26_east,data);
[Jac_cell_333,residual_cell_333] = problem_jac_res_init(q27,q27_east,data);



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


p = 
dpdn = 


Flux_x = @(q)[ q(2);
    q(2)^2/q(1) + p;
 (q(2)*q(3))/q(1);
 (q(2)*q(4))/q(1);
 (q(5)*q(2))/q(1);
 (q(6)*q(2))/q(1);
 (q(7)*q(2))/q(1);];
 
Flux_y = @(q)[ q(3);
 (q(2)*q(3))/q(1);
    q(3)^2/q(1) + p ;
 (q(3)*q(4))/q(1);
 (q(5)*q(3))/q(1);
 (q(6)*q(3))/q(1);
 (q(7)*q(3))/q(1);];
 
Flux_z = @(q)[ q(4);
 (q(2)*q(4))/q(1);
 (q(3)*q(4))/q(1);
    q(4)^2/q(1) + p ;
 (q(5)*q(4))/q(1);
 (q(6)*q(4))/q(1);
 (q(7)*q(4))/q(1);];


DFlux_dq_x = @(q) [[               0,        1,    0,    0,    0,    0,    0];
[ dpdn - q(2)^2/q(1)^2, (2*q(2))/q(1),    0,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2,     q(3)/q(1), q(2)/q(1),    0,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2,     q(4)/q(1),    0, q(2)/q(1),    0,    0,    0];
[    -(q(5)*q(2))/q(1)^2,     q(5)/q(1),    0,    0, q(2)/q(1),    0,    0];
[    -(q(6)*q(2))/q(1)^2,     q(6)/q(1),    0,    0,    0, q(2)/q(1),    0];
[    -(q(7)*q(2))/q(1)^2,     q(7)/q(1),    0,    0,    0,    0, q(2)/q(1)]];


DFlux_dq_y = @(q) [[               0,    0,        1,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2, q(3)/q(1),     q(2)/q(1),    0,    0,    0,    0];
[ dpdn - q(3)^2/q(1)^2,    0, (2*q(3))/q(1),    0,    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0,     q(4)/q(1), q(3)/q(1),    0,    0,    0];
[    -(q(5)*q(3))/q(1)^2,    0,     q(5)/q(1),    0, q(3)/q(1),    0,    0];
[    -(q(6)*q(3))/q(1)^2,    0,     q(6)/q(1),    0,    0, q(3)/q(1),    0];
[    -(q(7)*q(3))/q(1)^2,    0,     q(7)/q(1),    0,    0,    0, q(3)/q(1)]];

DFlux_dq_z = @(q) [[               0,    0,    0,        1,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2, q(4)/q(1),    0,     q(2)/q(1),    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0, q(4)/q(1),     q(3)/q(1),    0,    0,    0];
[ dpdn - q(4)^2/q(1)^2,    0,    0, (2*q(4))/q(1),    0,    0,    0];
[    -(q(5)*q(4))/q(1)^2,    0,    0,     q(5)/q(1), q(4)/q(1),    0,    0];
[    -(q(6)*q(4))/q(1)^2,    0,    0,     q(6)/q(1),    0, q(4)/q(1),    0];
[    -(q(7)*q(4))/q(1)^2,    0,    0,     q(7)/q(1),    0,    0, q(4)/q(1)];];


speed_x = @(q) abs(q(2)/q(1)) + abs(sqrt(dpdn)) ; 
speed_y = @(q) abs(q(3)/q(1)) + abs(sqrt(dpdn)) ; 
speed_z = @(q) abs(q(4)/q(1)) + abs(sqrt(dpdn)) ; 


for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = uadlocs(k,1);
    v1quad = Dlist(k,2);
    v2quad = Dlist(k,3);
    
    weight = uadwgts(k);    
    wgtnuv1 = weight*nuv1;
    wgtnuv2 = weight*nuv2;
    wgtnuv3 = weight*nuv3;
    
    psi_e = vectpsi_east(:,:,tquad,v1quad,v2quad);
    psi_w = vectpsi_west(:,:,tquad,v1quad,v2quad);
    psi_n = vectpsi_nort(:,:,tquad,v1quad,v2quad);
    psi_s = vectpsi_sout(:,:,tquad,v1quad,v2quad);
    psi_u = vectpsi_uppr(:,:,tquad,v1quad,v2quad);
    psi_d = vectpsi_down(:,:,tquad,v1quad,v2quad);
    
    psi_eT = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    psi_wT = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    psi_nT = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    psi_sT = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    psi_uT = vectpsi_uppr_Trans(:,:,tquad,v1quad,v2quad);
    psi_dT = vectpsi_down_Trans(:,:,tquad,v1quad,v2quad);
    
    % EAST AND WEST TERMS
    w01_e = psi_e*q01;
    w01_w = psi_w*q01;
    
    w02_e = psi_e*q02;
    w02_w = psi_w*q02;
    
    w03_e = psi_e*q03;
    w03_w = psi_w*q03;    
    
    w04_e = psi_e*q04;
    w04_w = psi_w*q04;    

    w05_e = psi_e*q05;
    w05_w = psi_w*q05;    

    w06_e = psi_e*q06;
    w06_w = psi_w*q06;    

    w07_e = psi_e*q07;
    w07_w = psi_w*q07;    

    w08_e = psi_e*q08;
    w08_w = psi_w*q08;    

    w09_e = psi_e*q09;
    w09_w = psi_w*q09;
    
    w10_e = psi_e*q10;
    w10_w = psi_w*q10;
    
    w11_e = psi_e*q11;
    w11_w = psi_w*q11;

    w12_e = psi_e*q12;
    w12_w = psi_w*q12;
   
    w13_e = psi_e*q13;
    w13_w = psi_w*q13;    
    
    w14_e = psi_e*q14;
    w14_w = psi_w*q14;    

    w15_e = psi_e*q15;
    w15_w = psi_w*q15;    

    w16_e = psi_e*q16;
    w16_w = psi_w*q16;    

    w17_e = psi_e*q17;
    w17_w = psi_w*q17;    

    w18_e = psi_e*q18;
    w18_w = psi_w*q18;    

    w19_e = psi_e*q19;
    w19_w = psi_w*q19;
    
    w20_e = psi_e*q20;
    w20_w = psi_w*q20;
    
    w21_e = psi_e*q21;
    w21_w = psi_w*q21;

    w22_e = psi_e*q22;
    w22_w = psi_w*q22;
   
    w23_e = psi_e*q23;
    w23_w = psi_w*q23;    
    
    w24_e = psi_e*q24;
    w24_w = psi_w*q24;    

    w25_e = psi_e*q25;
    w25_w = psi_w*q25;    

    w26_e = psi_e*q26;
    w26_w = psi_w*q26;    

    w27_e = psi_e*q27;
    w27_w = psi_w*q27;      

% NORT SOUT TERMS
    w01_n = psi_n*q01;
    w01_s = psi_s*q01;
    
    w02_n = psi_n*q02;
    w02_s = psi_s*q02;
    
    w03_n = psi_n*q03;
    w03_s = psi_s*q03;    
    
    w04_n = psi_n*q04;
    w04_s = psi_s*q04;    

    w05_n = psi_n*q05;
    w05_s = psi_s*q05;    

    w06_n = psi_n*q06;
    w06_s = psi_s*q06;    

    w07_n = psi_n*q07;
    w07_s = psi_s*q07;    

    w08_n = psi_n*q08;
    w08_s = psi_s*q08;    

    w09_n = psi_n*q09;
    w09_s = psi_s*q09;
    
    w10_n = psi_n*q10;
    w10_s = psi_s*q10;
    
    w11_n = psi_n*q11;
    w11_s = psi_s*q11;

    w12_n = psi_n*q12;
    w12_s = psi_s*q12;
   
    w13_n = psi_n*q13;
    w13_s = psi_s*q13;    
    
    w14_n = psi_n*q14;
    w14_s = psi_s*q14;    

    w15_n = psi_n*q15;
    w15_s = psi_s*q15;    

    w16_n = psi_n*q16;
    w16_s = psi_s*q16;    

    w17_n = psi_n*q17;
    w17_s = psi_s*q17;    

    w18_n = psi_n*q18;
    w18_s = psi_s*q18;    

    w19_n = psi_n*q19;
    w19_s = psi_s*q19;
    
    w20_n = psi_n*q20;
    w20_s = psi_s*q20;
    
    w21_n = psi_n*q21;
    w21_s = psi_s*q21;

    w22_n = psi_n*q22;
    w22_s = psi_s*q22;
   
    w23_n = psi_n*q23;
    w23_s = psi_s*q23;    
    
    w24_n = psi_n*q24;
    w24_s = psi_s*q24;    

    w25_n = psi_n*q25;
    w25_s = psi_s*q25;    

    w26_n = psi_n*q26;
    w26_s = psi_s*q26;    

    w27_n = psi_n*q27;
    w27_s = psi_s*q27;  
    
   

% UPPR DOWN TERMS
    w01_u = psi_u*q01;
    w01_d = psi_d*q01;
    
    w02_u = psi_u*q02;
    w02_d = psi_d*q02;
    
    w03_u = psi_u*q03;
    w03_d = psi_d*q03;    
    
    w04_u = psi_u*q04;
    w04_d = psi_d*q04;    

    w05_u = psi_u*q05;
    w05_d = psi_d*q05;    

    w06_u = psi_u*q06;
    w06_d = psi_d*q06;    

    w07_u = psi_u*q07;
    w07_d = psi_d*q07;    

    w08_u = psi_u*q08;
    w08_d = psi_d*q08;    

    w09_u = psi_u*q09;
    w09_d = psi_d*q09;
    
    w10_u = psi_u*q10;
    w10_d = psi_d*q10;
    
    w11_u = psi_u*q11;
    w11_d = psi_d*q11;

    w12_u = psi_u*q12;
    w12_d = psi_d*q12;
   
    w13_u = psi_u*q13;
    w13_d = psi_d*q13;    
    
    w14_u = psi_u*q14;
    w14_d = psi_d*q14;    

    w15_u = psi_u*q15;
    w15_d = psi_d*q15;    

    w16_u = psi_u*q16;
    w16_d = psi_d*q16;    

    w17_u = psi_u*q17;
    w17_d = psi_d*q17;    

    w18_u = psi_u*q18;
    w18_d = psi_d*q18;    

    w19_u = psi_u*q19;
    w19_d = psi_d*q19;
    
    w20_u = psi_u*q20;
    w20_d = psi_d*q20;
    
    w21_u = psi_u*q21;
    w21_d = psi_d*q21;

    w22_u = psi_u*q22;
    w22_d = psi_d*q22;
   
    w23_u = psi_u*q23;
    w23_d = psi_d*q23;    
    
    w24_u = psi_u*q24;
    w24_d = psi_d*q24;    

    w25_u = psi_u*q25;
    w25_d = psi_d*q25;    

    w26_u = psi_u*q26;
    w26_d = psi_d*q26;    

    w27_u = psi_u*q27;
    w27_d = psi_d*q27;  
        
    F1_e = Flux_x(w01_e);
    F2_e = Flux_x(w02_e);
    F3_e = Flux_x(w03_e);
    F4_e = Flux_x(w04_e);
    F5_e = Flux_x(w05_e);
    F6_e = Flux_x(w06_e);
    F7_e = Flux_x(w07_e);
    F8_e = Flux_x(w08_e);
    F9_e = Flux_x(w09_e);
    F10_e = Flux_x(w10_e);
    F11_e = Flux_x(w11_e);
    F12_e = Flux_x(w12_e);
    F13_e = Flux_x(w13_e);
    F14_e = Flux_x(w14_e);
    F15_e = Flux_x(w15_e);
    F16_e = Flux_x(w16_e);
    F17_e = Flux_x(w17_e);
    F18_e = Flux_x(w18_e);
    F19_e = Flux_x(w19_e);
    F20_e = Flux_x(w20_e);
    F21_e = Flux_x(w21_e);
    F22_e = Flux_x(w22_e);
    F23_e = Flux_x(w23_e);
    F24_e = Flux_x(w24_e);
    F25_e = Flux_x(w25_e);
    F26_e = Flux_x(w26_e);
    F27_e = Flux_x(w27_e);
    F1_w = Flux_x(w01_w);
    F2_w = Flux_x(w02_w);
    F3_w = Flux_x(w03_w);
    F4_w = Flux_x(w04_w);
    F5_w = Flux_x(w05_w);
    F6_w = Flux_x(w06_w);
    F7_w = Flux_x(w07_w);
    F8_w = Flux_x(w08_w);
    F9_w = Flux_x(w09_w);
    F10_w = Flux_x(w10_w);
    F11_w = Flux_x(w11_w);
    F12_w = Flux_x(w12_w);
    F13_w = Flux_x(w13_w);
    F14_w = Flux_x(w14_w);
    F15_w = Flux_x(w15_w);
    F16_w = Flux_x(w16_w);
    F17_w = Flux_x(w17_w);
    F18_w = Flux_x(w18_w);
    F19_w = Flux_x(w19_w);
    F20_w = Flux_x(w20_w);
    F21_w = Flux_x(w21_w);
    F22_w = Flux_x(w22_w);
    F23_w = Flux_x(w23_w);
    F24_w = Flux_x(w24_w);
    F25_w = Flux_x(w25_w);
    F26_w = Flux_x(w26_w);
    F27_w = Flux_x(w27_w);
    F1_n = Flux_y(w01_n);
    F2_n = Flux_y(w02_n);
    F3_n = Flux_y(w03_n);
    F4_n = Flux_y(w04_n);
    F5_n = Flux_y(w05_n);
    F6_n = Flux_y(w06_n);
    F7_n = Flux_y(w07_n);
    F8_n = Flux_y(w08_n);
    F9_n = Flux_y(w09_n);
    F10_n = Flux_y(w10_n);
    F11_n = Flux_y(w11_n);
    F12_n = Flux_y(w12_n);
    F13_n = Flux_y(w13_n);
    F14_n = Flux_y(w14_n);
    F15_n = Flux_y(w15_n);
    F16_n = Flux_y(w16_n);
    F17_n = Flux_y(w17_n);
    F18_n = Flux_y(w18_n);
    F19_n = Flux_y(w19_n);
    F20_n = Flux_y(w20_n);
    F21_n = Flux_y(w21_n);
    F22_n = Flux_y(w22_n);
    F23_n = Flux_y(w23_n);
    F24_n = Flux_y(w24_n);
    F25_n = Flux_y(w25_n);
    F26_n = Flux_y(w26_n);
    F27_n = Flux_y(w27_n);
    F1_s = Flux_y(w01_s);
    F2_s = Flux_y(w02_s);
    F3_s = Flux_y(w03_s);
    F4_s = Flux_y(w04_s);
    F5_s = Flux_y(w05_s);
    F6_s = Flux_y(w06_s);
    F7_s = Flux_y(w07_s);
    F8_s = Flux_y(w08_s);
    F9_s = Flux_y(w09_s);
    F10_s = Flux_y(w10_s);
    F11_s = Flux_y(w11_s);
    F12_s = Flux_y(w12_s);
    F13_s = Flux_y(w13_s);
    F14_s = Flux_y(w14_s);
    F15_s = Flux_y(w15_s);
    F16_s = Flux_y(w16_s);
    F17_s = Flux_y(w17_s);
    F18_s = Flux_y(w18_s);
    F19_s = Flux_y(w19_s);
    F20_s = Flux_y(w20_s);
    F21_s = Flux_y(w21_s);
    F22_s = Flux_y(w22_s);
    F23_s = Flux_y(w23_s);
    F24_s = Flux_y(w24_s);
    F25_s = Flux_y(w25_s);
    F26_s = Flux_y(w26_s);
    F27_s = Flux_y(w27_s);
    F1_u = Flux_z(w01_u);
    F2_u = Flux_z(w02_u);
    F3_u = Flux_z(w03_u);
    F4_u = Flux_z(w04_u);
    F5_u = Flux_z(w05_u);
    F6_u = Flux_z(w06_u);
    F7_u = Flux_z(w07_u);
    F8_u = Flux_z(w08_u);
    F9_u = Flux_z(w09_u);
    F10_u = Flux_z(w10_u);
    F11_u = Flux_z(w11_u);
    F12_u = Flux_z(w12_u);
    F13_u = Flux_z(w13_u);
    F14_u = Flux_z(w14_u);
    F15_u = Flux_z(w15_u);
    F16_u = Flux_z(w16_u);
    F17_u = Flux_z(w17_u);
    F18_u = Flux_z(w18_u);
    F19_u = Flux_z(w19_u);
    F20_u = Flux_z(w20_u);
    F21_u = Flux_z(w21_u);
    F22_u = Flux_z(w22_u);
    F23_u = Flux_z(w23_u);
    F24_u = Flux_z(w24_u);
    F25_u = Flux_z(w25_u);
    F26_u = Flux_z(w26_u);
    F27_u = Flux_z(w27_u);
    F1_d = Flux_z(w01_d);
    F2_d = Flux_z(w02_d);
    F3_d = Flux_z(w03_d);
    F4_d = Flux_z(w04_d);
    F5_d = Flux_z(w05_d);
    F6_d = Flux_z(w06_d);
    F7_d = Flux_z(w07_d);
    F8_d = Flux_z(w08_d);
    F9_d = Flux_z(w09_d);
    F10_d = Flux_z(w10_d);
    F11_d = Flux_z(w11_d);
    F12_d = Flux_z(w12_d);
    F13_d = Flux_z(w13_d);
    F14_d = Flux_z(w14_d);
    F15_d = Flux_z(w15_d);
    F16_d = Flux_z(w16_d);
    F17_d = Flux_z(w17_d);
    F18_d = Flux_z(w18_d);
    F19_d = Flux_z(w19_d);
    F20_d = Flux_z(w20_d);
    F21_d = Flux_z(w21_d);
    F22_d = Flux_z(w22_d);
    F23_d = Flux_z(w23_d);
    F24_d = Flux_z(w24_d);
    F25_d = Flux_z(w25_d);
    F26_d = Flux_z(w26_d);
    F27_d = Flux_z(w27_d);
         
    DF1_e = DFlux_dq_x(w01_e);
    DF2_e = DFlux_dq_x(w02_e);
    DF3_e = DFlux_dq_x(w03_e);
    DF4_e = DFlux_dq_x(w04_e);
    DF5_e = DFlux_dq_x(w05_e);
    DF6_e = DFlux_dq_x(w06_e);
    DF7_e = DFlux_dq_x(w07_e);
    DF8_e = DFlux_dq_x(w08_e);
    DF9_e = DFlux_dq_x(w09_e);
    DF10_e = DFlux_dq_x(w10_e);
    DF11_e = DFlux_dq_x(w11_e);
    DF12_e = DFlux_dq_x(w12_e);
    DF13_e = DFlux_dq_x(w13_e);
    DF14_e = DFlux_dq_x(w14_e);
    DF15_e = DFlux_dq_x(w15_e);
    DF16_e = DFlux_dq_x(w16_e);
    DF17_e = DFlux_dq_x(w17_e);
    DF18_e = DFlux_dq_x(w18_e);
    DF19_e = DFlux_dq_x(w19_e);
    DF20_e = DFlux_dq_x(w20_e);
    DF21_e = DFlux_dq_x(w21_e);
    DF22_e = DFlux_dq_x(w22_e);
    DF23_e = DFlux_dq_x(w23_e);
    DF24_e = DFlux_dq_x(w24_e);
    DF25_e = DFlux_dq_x(w25_e);
    DF26_e = DFlux_dq_x(w26_e);
    DF27_e = DFlux_dq_x(w27_e);
    DF1_w = DFlux_dq_x(w01_w);
    DF2_w = DFlux_dq_x(w02_w);
    DF3_w = DFlux_dq_x(w03_w);
    DF4_w = DFlux_dq_x(w04_w);
    DF5_w = DFlux_dq_x(w05_w);
    DF6_w = DFlux_dq_x(w06_w);
    DF7_w = DFlux_dq_x(w07_w);
    DF8_w = DFlux_dq_x(w08_w);
    DF9_w = DFlux_dq_x(w09_w);
    DF10_w = DFlux_dq_x(w10_w);
    DF11_w = DFlux_dq_x(w11_w);
    DF12_w = DFlux_dq_x(w12_w);
    DF13_w = DFlux_dq_x(w13_w);
    DF14_w = DFlux_dq_x(w14_w);
    DF15_w = DFlux_dq_x(w15_w);
    DF16_w = DFlux_dq_x(w16_w);
    DF17_w = DFlux_dq_x(w17_w);
    DF18_w = DFlux_dq_x(w18_w);
    DF19_w = DFlux_dq_x(w19_w);
    DF20_w = DFlux_dq_x(w20_w);
    DF21_w = DFlux_dq_x(w21_w);
    DF22_w = DFlux_dq_x(w22_w);
    DF23_w = DFlux_dq_x(w23_w);
    DF24_w = DFlux_dq_x(w24_w);
    DF25_w = DFlux_dq_x(w25_w);
    DF26_w = DFlux_dq_x(w26_w);
    DF27_w = DFlux_dq_x(w27_w);
    DF1_n = DFlux_dq_y(w01_n);
    DF2_n = DFlux_dq_y(w02_n);
    DF3_n = DFlux_dq_y(w03_n);
    DF4_n = DFlux_dq_y(w04_n);
    DF5_n = DFlux_dq_y(w05_n);
    DF6_n = DFlux_dq_y(w06_n);
    DF7_n = DFlux_dq_y(w07_n);
    DF8_n = DFlux_dq_y(w08_n);
    DF9_n = DFlux_dq_y(w09_n);
    DF10_n = DFlux_dq_y(w10_n);
    DF11_n = DFlux_dq_y(w11_n);
    DF12_n = DFlux_dq_y(w12_n);
    DF13_n = DFlux_dq_y(w13_n);
    DF14_n = DFlux_dq_y(w14_n);
    DF15_n = DFlux_dq_y(w15_n);
    DF16_n = DFlux_dq_y(w16_n);
    DF17_n = DFlux_dq_y(w17_n);
    DF18_n = DFlux_dq_y(w18_n);
    DF19_n = DFlux_dq_y(w19_n);
    DF20_n = DFlux_dq_y(w20_n);
    DF21_n = DFlux_dq_y(w21_n);
    DF22_n = DFlux_dq_y(w22_n);
    DF23_n = DFlux_dq_y(w23_n);
    DF24_n = DFlux_dq_y(w24_n);
    DF25_n = DFlux_dq_y(w25_n);
    DF26_n = DFlux_dq_y(w26_n);
    DF27_n = DFlux_dq_y(w27_n);
    DF1_s = DFlux_dq_y(w01_s);
    DF2_s = DFlux_dq_y(w02_s);
    DF3_s = DFlux_dq_y(w03_s);
    DF4_s = DFlux_dq_y(w04_s);
    DF5_s = DFlux_dq_y(w05_s);
    DF6_s = DFlux_dq_y(w06_s);
    DF7_s = DFlux_dq_y(w07_s);
    DF8_s = DFlux_dq_y(w08_s);
    DF9_s = DFlux_dq_y(w09_s);
    DF10_s = DFlux_dq_y(w10_s);
    DF11_s = DFlux_dq_y(w11_s);
    DF12_s = DFlux_dq_y(w12_s);
    DF13_s = DFlux_dq_y(w13_s);
    DF14_s = DFlux_dq_y(w14_s);
    DF15_s = DFlux_dq_y(w15_s);
    DF16_s = DFlux_dq_y(w16_s);
    DF17_s = DFlux_dq_y(w17_s);
    DF18_s = DFlux_dq_y(w18_s);
    DF19_s = DFlux_dq_y(w19_s);
    DF20_s = DFlux_dq_y(w20_s);
    DF21_s = DFlux_dq_y(w21_s);
    DF22_s = DFlux_dq_y(w22_s);
    DF23_s = DFlux_dq_y(w23_s);
    DF24_s = DFlux_dq_y(w24_s);
    DF25_s = DFlux_dq_y(w25_s);
    DF26_s = DFlux_dq_y(w26_s);
    DF27_s = DFlux_dq_y(w27_s);
    DF1_u = DFlux_dq_z(w01_u);
    DF2_u = DFlux_dq_z(w02_u);
    DF3_u = DFlux_dq_z(w03_u);
    DF4_u = DFlux_dq_z(w04_u);
    DF5_u = DFlux_dq_z(w05_u);
    DF6_u = DFlux_dq_z(w06_u);
    DF7_u = DFlux_dq_z(w07_u);
    DF8_u = DFlux_dq_z(w08_u);
    DF9_u = DFlux_dq_z(w09_u);
    DF10_u = DFlux_dq_z(w10_u);
    DF11_u = DFlux_dq_z(w11_u);
    DF12_u = DFlux_dq_z(w12_u);
    DF13_u = DFlux_dq_z(w13_u);
    DF14_u = DFlux_dq_z(w14_u);
    DF15_u = DFlux_dq_z(w15_u);
    DF16_u = DFlux_dq_z(w16_u);
    DF17_u = DFlux_dq_z(w17_u);
    DF18_u = DFlux_dq_z(w18_u);
    DF19_u = DFlux_dq_z(w19_u);
    DF20_u = DFlux_dq_z(w20_u);
    DF21_u = DFlux_dq_z(w21_u);
    DF22_u = DFlux_dq_z(w22_u);
    DF23_u = DFlux_dq_z(w23_u);
    DF24_u = DFlux_dq_z(w24_u);
    DF25_u = DFlux_dq_z(w25_u);
    DF26_u = DFlux_dq_z(w26_u);
    DF27_u = DFlux_dq_z(w27_u);
    DF1_d = DFlux_dq_z(w01_d);
    DF2_d = DFlux_dq_z(w02_d);
    DF3_d = DFlux_dq_z(w03_d);
    DF4_d = DFlux_dq_z(w04_d);
    DF5_d = DFlux_dq_z(w05_d);
    DF6_d = DFlux_dq_z(w06_d);
    DF7_d = DFlux_dq_z(w07_d);
    DF8_d = DFlux_dq_z(w08_d);
    DF9_d = DFlux_dq_z(w09_d);
    DF10_d = DFlux_dq_z(w10_d);
    DF11_d = DFlux_dq_z(w11_d);
    DF12_d = DFlux_dq_z(w12_d);
    DF13_d = DFlux_dq_z(w13_d);
    DF14_d = DFlux_dq_z(w14_d);
    DF15_d = DFlux_dq_z(w15_d);
    DF16_d = DFlux_dq_z(w16_d);
    DF17_d = DFlux_dq_z(w17_d);
    DF18_d = DFlux_dq_z(w18_d);
    DF19_d = DFlux_dq_z(w19_d);
    DF20_d = DFlux_dq_z(w20_d);
    DF21_d = DFlux_dq_z(w21_d);
    DF22_d = DFlux_dq_z(w22_d);
    DF23_d = DFlux_dq_z(w23_d);
    DF24_d = DFlux_dq_z(w24_d);
    DF25_d = DFlux_dq_z(w25_d);
    DF26_d = DFlux_dq_z(w26_d);
    DF27_d = DFlux_dq_z(w27_d);
    
    east_trunc_311 = east_trunc_311 + psi_eT*(F3_e*wgtnuv1);
    east_trunc_321 = east_trunc_321 + psi_eT*(F6_e*wgtnuv1);
    east_trunc_331 = east_trunc_331 + psi_eT*(F9_e*wgtnuv1);
    east_trunc_312 = east_trunc_312 + psi_eT*(F12_e*wgtnuv1);
    east_trunc_322 = east_trunc_322 + psi_eT*(F15_e*wgtnuv1);
    east_trunc_332 = east_trunc_332 + psi_eT*(F18_e*wgtnuv1);
    east_trunc_313 = east_trunc_313 + psi_eT*(F21_e*wgtnuv1);
    east_trunc_323 = east_trunc_323 + psi_eT*(F24_e*wgtnuv1);
    east_trunc_333 = east_trunc_333 + psi_eT*(F27_e*wgtnuv1);        
    
    Jac_east_trunc_311 = Jac_east_trunc_311 + psi_eT*DF3_e*psi_e*wgtnuv1;
    Jac_east_trunc_321 = Jac_east_trunc_321 + psi_eT*DF6_e*psi_e*wgtnuv1;
    Jac_east_trunc_331 = Jac_east_trunc_331 + psi_eT*DF9_e*psi_e*wgtnuv1;          
    Jac_east_trunc_312 = Jac_east_trunc_312 + psi_eT*DF12_e*psi_e*wgtnuv1;
    Jac_east_trunc_322 = Jac_east_trunc_322 + psi_eT*DF15_e*psi_e*wgtnuv1;
    Jac_east_trunc_332 = Jac_east_trunc_332 + psi_eT*DF18_e*psi_e*wgtnuv1;      
    Jac_east_trunc_313 = Jac_east_trunc_313 + psi_eT*DF21_e*psi_e*wgtnuv1;
    Jac_east_trunc_323 = Jac_east_trunc_323 + psi_eT*DF24_e*psi_e*wgtnuv1;
    Jac_east_trunc_333 = Jac_east_trunc_333 + psi_eT*DF27_e*psi_e*wgtnuv1;        
   
       
    west_trunc_111 = west_trunc_111 - psi_wT*(F1_w*wgtnuv1);
    west_trunc_121 = west_trunc_121 - psi_wT*(F4_w*wgtnuv1);
    west_trunc_131 = west_trunc_131 - psi_wT*(F7_w*wgtnuv1);           
    west_trunc_112 = west_trunc_112 - psi_wT*(F10_w*wgtnuv1);
    west_trunc_122 = west_trunc_122 - psi_wT*(F13_w*wgtnuv1);
    west_trunc_132 = west_trunc_132 - psi_wT*(F16_w*wgtnuv1);           
    west_trunc_113 = west_trunc_113 - psi_wT*(F19_w*wgtnuv1);
    west_trunc_123 = west_trunc_123 - psi_wT*(F22_w*wgtnuv1);
    west_trunc_133 = west_trunc_133 - psi_wT*(F25_w*wgtnuv1);
        
    Jac_west_trunc_111 = Jac_west_trunc_111 - psi_wT*(DF1_w*psi_w*wgtnuv1);
    Jac_west_trunc_121 = Jac_west_trunc_121 - psi_wT*(DF4_w*psi_w*wgtnuv1);
    Jac_west_trunc_131 = Jac_west_trunc_131 - psi_wT*(DF7_w*psi_w*wgtnuv1);         
    Jac_west_trunc_112 = Jac_west_trunc_112 - psi_wT*(DF10_w*psi_w*wgtnuv1);
    Jac_west_trunc_122 = Jac_west_trunc_122 - psi_wT*(DF13_w*psi_w*wgtnuv1);
    Jac_west_trunc_132 = Jac_west_trunc_132 - psi_wT*(DF16_w*psi_w*wgtnuv1);
    Jac_west_trunc_113 = Jac_west_trunc_113 - psi_wT*(DF19_w*psi_w*wgtnuv1);
    Jac_west_trunc_123 = Jac_west_trunc_123 - psi_wT*(DF22_w*psi_w*wgtnuv1);
    Jac_west_trunc_133 = Jac_west_trunc_133 - psi_wT*(DF25_w*psi_w*wgtnuv1);
    
    nort_trunc_131 = nort_trunc_131 + psi_nT*(F7_n*wgtnuv2);
    nort_trunc_231 = nort_trunc_231 + psi_nT*(F8_n*wgtnuv2);
    nort_trunc_331 = nort_trunc_331 + psi_nT*(F9_n*wgtnuv2);
    nort_trunc_132 = nort_trunc_132 + psi_nT*(F16_n*wgtnuv2);
    nort_trunc_232 = nort_trunc_232 + psi_nT*(F17_n*wgtnuv2);
    nort_trunc_332 = nort_trunc_332 + psi_nT*(F18_n*wgtnuv2);
    nort_trunc_133 = nort_trunc_133 + psi_nT*(F25_n*wgtnuv2);
    nort_trunc_233 = nort_trunc_233 + psi_nT*(F26_n*wgtnuv2);
    nort_trunc_333 = nort_trunc_333 + psi_nT*(F27_n*wgtnuv2);

    Jac_nort_trunc_131 = Jac_nort_trunc_131 + psi_nT*(DF7_n*psi_n*wgtnuv2);
    Jac_nort_trunc_231 = Jac_nort_trunc_231 + psi_nT*(DF8_n*psi_n*wgtnuv2);
    Jac_nort_trunc_331 = Jac_nort_trunc_331 + psi_nT*(DF9_n*psi_n*wgtnuv2);
    Jac_nort_trunc_132 = Jac_nort_trunc_132 + psi_nT*(DF16_n*psi_n*wgtnuv2);
    Jac_nort_trunc_232 = Jac_nort_trunc_232 + psi_nT*(DF17_n*psi_n*wgtnuv2);
    Jac_nort_trunc_332 = Jac_nort_trunc_332 + psi_nT*(DF18_n*psi_n*wgtnuv2);
    Jac_nort_trunc_133 = Jac_nort_trunc_133 + psi_nT*(DF25_n*psi_n*wgtnuv2);
    Jac_nort_trunc_233 = Jac_nort_trunc_233 + psi_nT*(DF26_n*psi_n*wgtnuv2);
    Jac_nort_trunc_333 = Jac_nort_trunc_333 + psi_nT*(DF27_n*psi_n*wgtnuv2);    
    
    sout_trunc_111 = sout_trunc_111 + psi_sT*(F1_s*wgtnuv2);
    sout_trunc_211 = sout_trunc_211 + psi_sT*(F2_s*wgtnuv2);
    sout_trunc_311 = sout_trunc_311 + psi_sT*(F3_s*wgtnuv2);
    sout_trunc_112 = sout_trunc_112 + psi_sT*(F10_s*wgtnuv2);
    sout_trunc_212 = sout_trunc_212 + psi_sT*(F11_s*wgtnuv2);
    sout_trunc_312 = sout_trunc_312 + psi_sT*(F12_s*wgtnuv2);
    sout_trunc_113 = sout_trunc_113 + psi_sT*(F19_s*wgtnuv2);
    sout_trunc_213 = sout_trunc_213 + psi_sT*(F20_s*wgtnuv2);
    sout_trunc_313 = sout_trunc_313 + psi_sT*(F21_s*wgtnuv2);

    Jac_sout_trunc_111 = Jac_sout_trunc_111 + psi_sT*(DF1_s*psi_s*wgtnuv2);
    Jac_sout_trunc_211 = Jac_sout_trunc_211 + psi_sT*(DF2_s*psi_s*wgtnuv2);
    Jac_sout_trunc_311 = Jac_sout_trunc_311 + psi_sT*(DF3_s*psi_s*wgtnuv2);
    Jac_sout_trunc_112 = Jac_sout_trunc_112 + psi_sT*(DF10_s*psi_s*wgtnuv2);
    Jac_sout_trunc_212 = Jac_sout_trunc_212 + psi_sT*(DF11_s*psi_s*wgtnuv2);
    Jac_sout_trunc_312 = Jac_sout_trunc_312 + psi_sT*(DF12_s*psi_s*wgtnuv2);
    Jac_sout_trunc_113 = Jac_sout_trunc_113 + psi_sT*(DF19_s*psi_s*wgtnuv2);
    Jac_sout_trunc_213 = Jac_sout_trunc_213 + psi_sT*(DF20_s*psi_s*wgtnuv2);
    Jac_sout_trunc_313 = Jac_sout_trunc_313 + psi_sT*(DF21_s*psi_s*wgtnuv2);    
        
    uppr_trunc_213  = uppr_trunc_213  + psi_uT*(F20_u*wgtnuv3);
    uppr_trunc_313  = uppr_trunc_313  + psi_uT*(F21_u*wgtnuv3);
    uppr_trunc_123  = uppr_trunc_123  + psi_uT*(F22_u*wgtnuv3);
    uppr_trunc_223  = uppr_trunc_223  + psi_uT*(F23_u*wgtnuv3);
    uppr_trunc_323  = uppr_trunc_323  + psi_uT*(F24_u*wgtnuv3);
    uppr_trunc_133  = uppr_trunc_133  + psi_uT*(F25_u*wgtnuv3);
    uppr_trunc_233  = uppr_trunc_233  + psi_uT*(F26_u*wgtnuv3);
    uppr_trunc_333  = uppr_trunc_333  + psi_uT*(F27_u*wgtnuv3);

    Jac_uppr_trunc_113 = Jac_uppr_trunc_113 + psi_uT*(DF19_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_213 = Jac_uppr_trunc_213 + psi_uT*(DF20_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_313 = Jac_uppr_trunc_313 + psi_uT*(DF21_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_123 = Jac_uppr_trunc_123 + psi_uT*(DF22_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_223 = Jac_uppr_trunc_223 + psi_uT*(DF23_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_323 = Jac_uppr_trunc_323 + psi_uT*(DF24_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_133 = Jac_uppr_trunc_133 + psi_uT*(DF25_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_233 = Jac_uppr_trunc_233 + psi_uT*(DF26_u*psi_u*wgtnuv3);
    Jac_uppr_trunc_333 = Jac_uppr_trunc_333 + psi_uT*(DF27_u*psi_u*wgtnuv3);


    down_trunc_111  = down_trunc_111  + psi_dT*(F1_d*wgtnuv3);
    down_trunc_211  = down_trunc_211  + psi_dT*(F2_d*wgtnuv3);
    down_trunc_311  = down_trunc_311  + psi_dT*(F3_d*wgtnuv3);
    down_trunc_121  = down_trunc_121  + psi_dT*(F4_d*wgtnuv3);
    down_trunc_221  = down_trunc_221  + psi_dT*(F5_d*wgtnuv3);
    down_trunc_321  = down_trunc_321  + psi_dT*(F6_d*wgtnuv3);
    down_trunc_131  = down_trunc_131  + psi_dT*(F7_d*wgtnuv3);
    down_trunc_231  = down_trunc_231  + psi_dT*(F8_d*wgtnuv3);
    down_trunc_331  = down_trunc_331  + psi_dT*(F9_d*wgtnuv3);

    Jac_down_trunc_111 = Jac_down_trunc_111 + psi_dT*(DF1_d*psi_d*wgtnuv3);
    Jac_down_trunc_211 = Jac_down_trunc_211 + psi_dT*(DF2_d*psi_d*wgtnuv3);
    Jac_down_trunc_311 = Jac_down_trunc_311 + psi_dT*(DF3_d*psi_d*wgtnuv3);
    Jac_down_trunc_121 = Jac_down_trunc_121 + psi_dT*(DF4_d*psi_d*wgtnuv3);
    Jac_down_trunc_221 = Jac_down_trunc_221 + psi_dT*(DF5_d*psi_d*wgtnuv3);
    Jac_down_trunc_321 = Jac_down_trunc_321 + psi_dT*(DF6_d*psi_d*wgtnuv3);
    Jac_down_trunc_131 = Jac_down_trunc_131 + psi_dT*(DF7_d*psi_d*wgtnuv3);
    Jac_down_trunc_231 = Jac_down_trunc_231 + psi_dT*(DF8_d*psi_d*wgtnuv3);
    Jac_down_trunc_331 = Jac_down_trunc_331 + psi_dT*(DF9_d*psi_d*wgtnuv3);
    
    %Interface1-2
    %cellleft111
    %cellrite211
    %e1w
    %eastwest
    speedleft=speed_x(w1_e);
    speedavg=speed_x((w1_e+w2_w)/2);
    speedrite=speed_x(w2_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F1_e+F2_w-speed*(w2_w-w1_e))*wgtnuv1/2;
    F_dql=(DF1_e+speed)*wgtnuv1/2;
    F_dqr=(DF2_w-speed)*wgtnuv1/2;
    residual_cell_111=residual_cell_111+psieT*Flux;
    residual_cell_211=residual_cell_211-psiwT*Flux;
    Jac_cell_111=Jac_cell_111+psieT*F_dql*psie;
    Jac_east_othr_111=Jac_east_othr_111+psieT*F_dqr*psiw;
    Jac_cell_211=Jac_cell_211-psiwT*F_dqr*psiw;
    Jac_west_othr_211=Jac_west_othr_211-psiwT*F_dql*psie;

    %Interface2-3
    %cellleft211
    %cellrite311
    %e1w
    %eastwest
    speedleft=speed_x(w2_e);
    speedavg=speed_x((w2_e+w3_w)/2);
    speedrite=speed_x(w3_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F2_e+F3_w-speed*(w3_w-w2_e))*wgtnuv1/2;
    F_dql=(DF2_e+speed)*wgtnuv1/2;
    F_dqr=(DF3_w-speed)*wgtnuv1/2;
    residual_cell_211=residual_cell_211+psieT*Flux;
    residual_cell_311=residual_cell_311-psiwT*Flux;
    Jac_cell_211=Jac_cell_211+psieT*F_dql*psie;
    Jac_east_othr_211=Jac_east_othr_211+psieT*F_dqr*psiw;
    Jac_cell_311=Jac_cell_311-psiwT*F_dqr*psiw;
    Jac_west_othr_311=Jac_west_othr_311-psiwT*F_dql*psie;

    %Interface4-5
    %cellleft121
    %cellrite221
    %e1w
    %eastwest
    speedleft=speed_x(w4_e);
    speedavg=speed_x((w4_e+w5_w)/2);
    speedrite=speed_x(w5_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F4_e+F5_w-speed*(w5_w-w4_e))*wgtnuv1/2;
    F_dql=(DF4_e+speed)*wgtnuv1/2;
    F_dqr=(DF5_w-speed)*wgtnuv1/2;
    residual_cell_121=residual_cell_121+psieT*Flux;
    residual_cell_221=residual_cell_221-psiwT*Flux;
    Jac_cell_121=Jac_cell_121+psieT*F_dql*psie;
    Jac_east_othr_121=Jac_east_othr_121+psieT*F_dqr*psiw;
    Jac_cell_221=Jac_cell_221-psiwT*F_dqr*psiw;
    Jac_west_othr_221=Jac_west_othr_221-psiwT*F_dql*psie;

    %Interface5-6
    %cellleft221
    %cellrite321
    %e1w
    %eastwest
    speedleft=speed_x(w5_e);
    speedavg=speed_x((w5_e+w6_w)/2);
    speedrite=speed_x(w6_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F5_e+F6_w-speed*(w6_w-w5_e))*wgtnuv1/2;
    F_dql=(DF5_e+speed)*wgtnuv1/2;
    F_dqr=(DF6_w-speed)*wgtnuv1/2;
    residual_cell_221=residual_cell_221+psieT*Flux;
    residual_cell_321=residual_cell_321-psiwT*Flux;
    Jac_cell_221=Jac_cell_221+psieT*F_dql*psie;
    Jac_east_othr_221=Jac_east_othr_221+psieT*F_dqr*psiw;
    Jac_cell_321=Jac_cell_321-psiwT*F_dqr*psiw;
    Jac_west_othr_321=Jac_west_othr_321-psiwT*F_dql*psie;

    %Interface7-8
    %cellleft131
    %cellrite231
    %e1w
    %eastwest
    speedleft=speed_x(w7_e);
    speedavg=speed_x((w7_e+w8_w)/2);
    speedrite=speed_x(w8_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F7_e+F8_w-speed*(w8_w-w7_e))*wgtnuv1/2;
    F_dql=(DF7_e+speed)*wgtnuv1/2;
    F_dqr=(DF8_w-speed)*wgtnuv1/2;
    residual_cell_131=residual_cell_131+psieT*Flux;
    residual_cell_231=residual_cell_231-psiwT*Flux;
    Jac_cell_131=Jac_cell_131+psieT*F_dql*psie;
    Jac_east_othr_131=Jac_east_othr_131+psieT*F_dqr*psiw;
    Jac_cell_231=Jac_cell_231-psiwT*F_dqr*psiw;
    Jac_west_othr_231=Jac_west_othr_231-psiwT*F_dql*psie;

    %Interface8-9
    %cellleft231
    %cellrite331
    %e1w
    %eastwest
    speedleft=speed_x(w8_e);
    speedavg=speed_x((w8_e+w9_w)/2);
    speedrite=speed_x(w9_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F8_e+F9_w-speed*(w9_w-w8_e))*wgtnuv1/2;
    F_dql=(DF8_e+speed)*wgtnuv1/2;
    F_dqr=(DF9_w-speed)*wgtnuv1/2;
    residual_cell_231=residual_cell_231+psieT*Flux;
    residual_cell_331=residual_cell_331-psiwT*Flux;
    Jac_cell_231=Jac_cell_231+psieT*F_dql*psie;
    Jac_east_othr_231=Jac_east_othr_231+psieT*F_dqr*psiw;
    Jac_cell_331=Jac_cell_331-psiwT*F_dqr*psiw;
    Jac_west_othr_331=Jac_west_othr_331-psiwT*F_dql*psie;

    %Interface10-11
    %cellleft112
    %cellrite212
    %e1w
    %eastwest
    speedleft=speed_x(w10_e);
    speedavg=speed_x((w10_e+w11_w)/2);
    speedrite=speed_x(w11_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F10_e+F11_w-speed*(w11_w-w10_e))*wgtnuv1/2;
    F_dql=(DF10_e+speed)*wgtnuv1/2;
    F_dqr=(DF11_w-speed)*wgtnuv1/2;
    residual_cell_112=residual_cell_112+psieT*Flux;
    residual_cell_212=residual_cell_212-psiwT*Flux;
    Jac_cell_112=Jac_cell_112+psieT*F_dql*psie;
    Jac_east_othr_112=Jac_east_othr_112+psieT*F_dqr*psiw;
    Jac_cell_212=Jac_cell_212-psiwT*F_dqr*psiw;
    Jac_west_othr_212=Jac_west_othr_212-psiwT*F_dql*psie;

    %Interface11-12
    %cellleft212
    %cellrite312
    %e1w
    %eastwest
    speedleft=speed_x(w11_e);
    speedavg=speed_x((w11_e+w12_w)/2);
    speedrite=speed_x(w12_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F11_e+F12_w-speed*(w12_w-w11_e))*wgtnuv1/2;
    F_dql=(DF11_e+speed)*wgtnuv1/2;
    F_dqr=(DF12_w-speed)*wgtnuv1/2;
    residual_cell_212=residual_cell_212+psieT*Flux;
    residual_cell_312=residual_cell_312-psiwT*Flux;
    Jac_cell_212=Jac_cell_212+psieT*F_dql*psie;
    Jac_east_othr_212=Jac_east_othr_212+psieT*F_dqr*psiw;
    Jac_cell_312=Jac_cell_312-psiwT*F_dqr*psiw;
    Jac_west_othr_312=Jac_west_othr_312-psiwT*F_dql*psie;

    %Interface13-14
    %cellleft122
    %cellrite222
    %e1w
    %eastwest
    speedleft=speed_x(w13_e);
    speedavg=speed_x((w13_e+w14_w)/2);
    speedrite=speed_x(w14_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F13_e+F14_w-speed*(w14_w-w13_e))*wgtnuv1/2;
    F_dql=(DF13_e+speed)*wgtnuv1/2;
    F_dqr=(DF14_w-speed)*wgtnuv1/2;
    residual_cell_122=residual_cell_122+psieT*Flux;
    residual_cell_222=residual_cell_222-psiwT*Flux;
    Jac_cell_122=Jac_cell_122+psieT*F_dql*psie;
    Jac_east_othr_122=Jac_east_othr_122+psieT*F_dqr*psiw;
    Jac_cell_222=Jac_cell_222-psiwT*F_dqr*psiw;
    Jac_west_othr_222=Jac_west_othr_222-psiwT*F_dql*psie;

    %Interface14-15
    %cellleft222
    %cellrite322
    %e1w
    %eastwest
    speedleft=speed_x(w14_e);
    speedavg=speed_x((w14_e+w15_w)/2);
    speedrite=speed_x(w15_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F14_e+F15_w-speed*(w15_w-w14_e))*wgtnuv1/2;
    F_dql=(DF14_e+speed)*wgtnuv1/2;
    F_dqr=(DF15_w-speed)*wgtnuv1/2;
    residual_cell_222=residual_cell_222+psieT*Flux;
    residual_cell_322=residual_cell_322-psiwT*Flux;
    Jac_cell_222=Jac_cell_222+psieT*F_dql*psie;
    Jac_east_othr_222=Jac_east_othr_222+psieT*F_dqr*psiw;
    Jac_cell_322=Jac_cell_322-psiwT*F_dqr*psiw;
    Jac_west_othr_322=Jac_west_othr_322-psiwT*F_dql*psie;

    %Interface16-17
    %cellleft132
    %cellrite232
    %e1w
    %eastwest
    speedleft=speed_x(w16_e);
    speedavg=speed_x((w16_e+w17_w)/2);
    speedrite=speed_x(w17_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F16_e+F17_w-speed*(w17_w-w16_e))*wgtnuv1/2;
    F_dql=(DF16_e+speed)*wgtnuv1/2;
    F_dqr=(DF17_w-speed)*wgtnuv1/2;
    residual_cell_132=residual_cell_132+psieT*Flux;
    residual_cell_232=residual_cell_232-psiwT*Flux;
    Jac_cell_132=Jac_cell_132+psieT*F_dql*psie;
    Jac_east_othr_132=Jac_east_othr_132+psieT*F_dqr*psiw;
    Jac_cell_232=Jac_cell_232-psiwT*F_dqr*psiw;
    Jac_west_othr_232=Jac_west_othr_232-psiwT*F_dql*psie;

    %Interface17-18
    %cellleft232
    %cellrite332
    %e1w
    %eastwest
    speedleft=speed_x(w17_e);
    speedavg=speed_x((w17_e+w18_w)/2);
    speedrite=speed_x(w18_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F17_e+F18_w-speed*(w18_w-w17_e))*wgtnuv1/2;
    F_dql=(DF17_e+speed)*wgtnuv1/2;
    F_dqr=(DF18_w-speed)*wgtnuv1/2;
    residual_cell_232=residual_cell_232+psieT*Flux;
    residual_cell_332=residual_cell_332-psiwT*Flux;
    Jac_cell_232=Jac_cell_232+psieT*F_dql*psie;
    Jac_east_othr_232=Jac_east_othr_232+psieT*F_dqr*psiw;
    Jac_cell_332=Jac_cell_332-psiwT*F_dqr*psiw;
    Jac_west_othr_332=Jac_west_othr_332-psiwT*F_dql*psie;

    %Interface19-20
    %cellleft113
    %cellrite213
    %e1w
    %eastwest
    speedleft=speed_x(w19_e);
    speedavg=speed_x((w19_e+w20_w)/2);
    speedrite=speed_x(w20_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F19_e+F20_w-speed*(w20_w-w19_e))*wgtnuv1/2;
    F_dql=(DF19_e+speed)*wgtnuv1/2;
    F_dqr=(DF20_w-speed)*wgtnuv1/2;
    residual_cell_113=residual_cell_113+psieT*Flux;
    residual_cell_213=residual_cell_213-psiwT*Flux;
    Jac_cell_113=Jac_cell_113+psieT*F_dql*psie;
    Jac_east_othr_113=Jac_east_othr_113+psieT*F_dqr*psiw;
    Jac_cell_213=Jac_cell_213-psiwT*F_dqr*psiw;
    Jac_west_othr_213=Jac_west_othr_213-psiwT*F_dql*psie;

    %Interface20-21
    %cellleft213
    %cellrite313
    %e1w
    %eastwest
    speedleft=speed_x(w20_e);
    speedavg=speed_x((w20_e+w21_w)/2);
    speedrite=speed_x(w21_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F20_e+F21_w-speed*(w21_w-w20_e))*wgtnuv1/2;
    F_dql=(DF20_e+speed)*wgtnuv1/2;
    F_dqr=(DF21_w-speed)*wgtnuv1/2;
    residual_cell_213=residual_cell_213+psieT*Flux;
    residual_cell_313=residual_cell_313-psiwT*Flux;
    Jac_cell_213=Jac_cell_213+psieT*F_dql*psie;
    Jac_east_othr_213=Jac_east_othr_213+psieT*F_dqr*psiw;
    Jac_cell_313=Jac_cell_313-psiwT*F_dqr*psiw;
    Jac_west_othr_313=Jac_west_othr_313-psiwT*F_dql*psie;

    %Interface22-23
    %cellleft123
    %cellrite223
    %e1w
    %eastwest
    speedleft=speed_x(w22_e);
    speedavg=speed_x((w22_e+w23_w)/2);
    speedrite=speed_x(w23_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F22_e+F23_w-speed*(w23_w-w22_e))*wgtnuv1/2;
    F_dql=(DF22_e+speed)*wgtnuv1/2;
    F_dqr=(DF23_w-speed)*wgtnuv1/2;
    residual_cell_123=residual_cell_123+psieT*Flux;
    residual_cell_223=residual_cell_223-psiwT*Flux;
    Jac_cell_123=Jac_cell_123+psieT*F_dql*psie;
    Jac_east_othr_123=Jac_east_othr_123+psieT*F_dqr*psiw;
    Jac_cell_223=Jac_cell_223-psiwT*F_dqr*psiw;
    Jac_west_othr_223=Jac_west_othr_223-psiwT*F_dql*psie;

    %Interface23-24
    %cellleft223
    %cellrite323
    %e1w
    %eastwest
    speedleft=speed_x(w23_e);
    speedavg=speed_x((w23_e+w24_w)/2);
    speedrite=speed_x(w24_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F23_e+F24_w-speed*(w24_w-w23_e))*wgtnuv1/2;
    F_dql=(DF23_e+speed)*wgtnuv1/2;
    F_dqr=(DF24_w-speed)*wgtnuv1/2;
    residual_cell_223=residual_cell_223+psieT*Flux;
    residual_cell_323=residual_cell_323-psiwT*Flux;
    Jac_cell_223=Jac_cell_223+psieT*F_dql*psie;
    Jac_east_othr_223=Jac_east_othr_223+psieT*F_dqr*psiw;
    Jac_cell_323=Jac_cell_323-psiwT*F_dqr*psiw;
    Jac_west_othr_323=Jac_west_othr_323-psiwT*F_dql*psie;

    %Interface25-26
    %cellleft133
    %cellrite233
    %e1w
    %eastwest
    speedleft=speed_x(w25_e);
    speedavg=speed_x((w25_e+w26_w)/2);
    speedrite=speed_x(w26_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F25_e+F26_w-speed*(w26_w-w25_e))*wgtnuv1/2;
    F_dql=(DF25_e+speed)*wgtnuv1/2;
    F_dqr=(DF26_w-speed)*wgtnuv1/2;
    residual_cell_133=residual_cell_133+psieT*Flux;
    residual_cell_233=residual_cell_233-psiwT*Flux;
    Jac_cell_133=Jac_cell_133+psieT*F_dql*psie;
    Jac_east_othr_133=Jac_east_othr_133+psieT*F_dqr*psiw;
    Jac_cell_233=Jac_cell_233-psiwT*F_dqr*psiw;
    Jac_west_othr_233=Jac_west_othr_233-psiwT*F_dql*psie;

    %Interface26-27
    %cellleft233
    %cellrite333
    %e1w
    %eastwest
    speedleft=speed_x(w26_e);
    speedavg=speed_x((w26_e+w27_w)/2);
    speedrite=speed_x(w27_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F26_e+F27_w-speed*(w27_w-w26_e))*wgtnuv1/2;
    F_dql=(DF26_e+speed)*wgtnuv1/2;
    F_dqr=(DF27_w-speed)*wgtnuv1/2;
    residual_cell_233=residual_cell_233+psieT*Flux;
    residual_cell_333=residual_cell_333-psiwT*Flux;
    Jac_cell_233=Jac_cell_233+psieT*F_dql*psie;
    Jac_east_othr_233=Jac_east_othr_233+psieT*F_dqr*psiw;
    Jac_cell_333=Jac_cell_333-psiwT*F_dqr*psiw;
    Jac_west_othr_333=Jac_west_othr_333-psiwT*F_dql*psie;


    %Interface1-4
    %cellleft111
    %cellrite121
    %n2s
    %nortsout
    speedleft=speed_y(w1_n);
    speedavg=speed_y((w1_n+w4_s)/2);
    speedrite=speed_y(w4_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F1_n+F4_s-speed*(w4_s-w1_n))*wgtnuv2/2;
    F_dql=(DF1_n+speed)*wgtnuv2/2;
    F_dqr=(DF4_s-speed)*wgtnuv2/2;
    residual_cell_111=residual_cell_111+psinT*Flux;
    residual_cell_121=residual_cell_121-psisT*Flux;
    Jac_cell_111=Jac_cell_111+psinT*F_dql*psin;
    Jac_nort_othr_111=Jac_nort_othr_111+psinT*F_dqr*psis;
    Jac_cell_121=Jac_cell_121-psisT*F_dqr*psis;
    Jac_sout_othr_121=Jac_sout_othr_121-psisT*F_dql*psin;

    %Interface4-7
    %cellleft121
    %cellrite131
    %n2s
    %nortsout
    speedleft=speed_y(w4_n);
    speedavg=speed_y((w4_n+w7_s)/2);
    speedrite=speed_y(w7_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F4_n+F7_s-speed*(w7_s-w4_n))*wgtnuv2/2;
    F_dql=(DF4_n+speed)*wgtnuv2/2;
    F_dqr=(DF7_s-speed)*wgtnuv2/2;
    residual_cell_121=residual_cell_121+psinT*Flux;
    residual_cell_131=residual_cell_131-psisT*Flux;
    Jac_cell_121=Jac_cell_121+psinT*F_dql*psin;
    Jac_nort_othr_121=Jac_nort_othr_121+psinT*F_dqr*psis;
    Jac_cell_131=Jac_cell_131-psisT*F_dqr*psis;
    Jac_sout_othr_131=Jac_sout_othr_131-psisT*F_dql*psin;

    %Interface2-5
    %cellleft211
    %cellrite221
    %n2s
    %nortsout
    speedleft=speed_y(w2_n);
    speedavg=speed_y((w2_n+w5_s)/2);
    speedrite=speed_y(w5_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F2_n+F5_s-speed*(w5_s-w2_n))*wgtnuv2/2;
    F_dql=(DF2_n+speed)*wgtnuv2/2;
    F_dqr=(DF5_s-speed)*wgtnuv2/2;
    residual_cell_211=residual_cell_211+psinT*Flux;
    residual_cell_221=residual_cell_221-psisT*Flux;
    Jac_cell_211=Jac_cell_211+psinT*F_dql*psin;
    Jac_nort_othr_211=Jac_nort_othr_211+psinT*F_dqr*psis;
    Jac_cell_221=Jac_cell_221-psisT*F_dqr*psis;
    Jac_sout_othr_221=Jac_sout_othr_221-psisT*F_dql*psin;

    %Interface5-8
    %cellleft221
    %cellrite231
    %n2s
    %nortsout
    speedleft=speed_y(w5_n);
    speedavg=speed_y((w5_n+w8_s)/2);
    speedrite=speed_y(w8_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F5_n+F8_s-speed*(w8_s-w5_n))*wgtnuv2/2;
    F_dql=(DF5_n+speed)*wgtnuv2/2;
    F_dqr=(DF8_s-speed)*wgtnuv2/2;
    residual_cell_221=residual_cell_221+psinT*Flux;
    residual_cell_231=residual_cell_231-psisT*Flux;
    Jac_cell_221=Jac_cell_221+psinT*F_dql*psin;
    Jac_nort_othr_221=Jac_nort_othr_221+psinT*F_dqr*psis;
    Jac_cell_231=Jac_cell_231-psisT*F_dqr*psis;
    Jac_sout_othr_231=Jac_sout_othr_231-psisT*F_dql*psin;

    %Interface3-6
    %cellleft311
    %cellrite321
    %n2s
    %nortsout
    speedleft=speed_y(w3_n);
    speedavg=speed_y((w3_n+w6_s)/2);
    speedrite=speed_y(w6_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F3_n+F6_s-speed*(w6_s-w3_n))*wgtnuv2/2;
    F_dql=(DF3_n+speed)*wgtnuv2/2;
    F_dqr=(DF6_s-speed)*wgtnuv2/2;
    residual_cell_311=residual_cell_311+psinT*Flux;
    residual_cell_321=residual_cell_321-psisT*Flux;
    Jac_cell_311=Jac_cell_311+psinT*F_dql*psin;
    Jac_nort_othr_311=Jac_nort_othr_311+psinT*F_dqr*psis;
    Jac_cell_321=Jac_cell_321-psisT*F_dqr*psis;
    Jac_sout_othr_321=Jac_sout_othr_321-psisT*F_dql*psin;

    %Interface6-9
    %cellleft321
    %cellrite331
    %n2s
    %nortsout
    speedleft=speed_y(w6_n);
    speedavg=speed_y((w6_n+w9_s)/2);
    speedrite=speed_y(w9_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F6_n+F9_s-speed*(w9_s-w6_n))*wgtnuv2/2;
    F_dql=(DF6_n+speed)*wgtnuv2/2;
    F_dqr=(DF9_s-speed)*wgtnuv2/2;
    residual_cell_321=residual_cell_321+psinT*Flux;
    residual_cell_331=residual_cell_331-psisT*Flux;
    Jac_cell_321=Jac_cell_321+psinT*F_dql*psin;
    Jac_nort_othr_321=Jac_nort_othr_321+psinT*F_dqr*psis;
    Jac_cell_331=Jac_cell_331-psisT*F_dqr*psis;
    Jac_sout_othr_331=Jac_sout_othr_331-psisT*F_dql*psin;

    %Interface10-13
    %cellleft112
    %cellrite122
    %n2s
    %nortsout
    speedleft=speed_y(w10_n);
    speedavg=speed_y((w10_n+w13_s)/2);
    speedrite=speed_y(w13_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F10_n+F13_s-speed*(w13_s-w10_n))*wgtnuv2/2;
    F_dql=(DF10_n+speed)*wgtnuv2/2;
    F_dqr=(DF13_s-speed)*wgtnuv2/2;
    residual_cell_112=residual_cell_112+psinT*Flux;
    residual_cell_122=residual_cell_122-psisT*Flux;
    Jac_cell_112=Jac_cell_112+psinT*F_dql*psin;
    Jac_nort_othr_112=Jac_nort_othr_112+psinT*F_dqr*psis;
    Jac_cell_122=Jac_cell_122-psisT*F_dqr*psis;
    Jac_sout_othr_122=Jac_sout_othr_122-psisT*F_dql*psin;

    %Interface13-16
    %cellleft122
    %cellrite132
    %n2s
    %nortsout
    speedleft=speed_y(w13_n);
    speedavg=speed_y((w13_n+w16_s)/2);
    speedrite=speed_y(w16_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F13_n+F16_s-speed*(w16_s-w13_n))*wgtnuv2/2;
    F_dql=(DF13_n+speed)*wgtnuv2/2;
    F_dqr=(DF16_s-speed)*wgtnuv2/2;
    residual_cell_122=residual_cell_122+psinT*Flux;
    residual_cell_132=residual_cell_132-psisT*Flux;
    Jac_cell_122=Jac_cell_122+psinT*F_dql*psin;
    Jac_nort_othr_122=Jac_nort_othr_122+psinT*F_dqr*psis;
    Jac_cell_132=Jac_cell_132-psisT*F_dqr*psis;
    Jac_sout_othr_132=Jac_sout_othr_132-psisT*F_dql*psin;

    %Interface11-14
    %cellleft212
    %cellrite222
    %n2s
    %nortsout
    speedleft=speed_y(w11_n);
    speedavg=speed_y((w11_n+w14_s)/2);
    speedrite=speed_y(w14_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F11_n+F14_s-speed*(w14_s-w11_n))*wgtnuv2/2;
    F_dql=(DF11_n+speed)*wgtnuv2/2;
    F_dqr=(DF14_s-speed)*wgtnuv2/2;
    residual_cell_212=residual_cell_212+psinT*Flux;
    residual_cell_222=residual_cell_222-psisT*Flux;
    Jac_cell_212=Jac_cell_212+psinT*F_dql*psin;
    Jac_nort_othr_212=Jac_nort_othr_212+psinT*F_dqr*psis;
    Jac_cell_222=Jac_cell_222-psisT*F_dqr*psis;
    Jac_sout_othr_222=Jac_sout_othr_222-psisT*F_dql*psin;

    %Interface14-17
    %cellleft222
    %cellrite232
    %n2s
    %nortsout
    speedleft=speed_y(w14_n);
    speedavg=speed_y((w14_n+w17_s)/2);
    speedrite=speed_y(w17_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F14_n+F17_s-speed*(w17_s-w14_n))*wgtnuv2/2;
    F_dql=(DF14_n+speed)*wgtnuv2/2;
    F_dqr=(DF17_s-speed)*wgtnuv2/2;
    residual_cell_222=residual_cell_222+psinT*Flux;
    residual_cell_232=residual_cell_232-psisT*Flux;
    Jac_cell_222=Jac_cell_222+psinT*F_dql*psin;
    Jac_nort_othr_222=Jac_nort_othr_222+psinT*F_dqr*psis;
    Jac_cell_232=Jac_cell_232-psisT*F_dqr*psis;
    Jac_sout_othr_232=Jac_sout_othr_232-psisT*F_dql*psin;

    %Interface12-15
    %cellleft312
    %cellrite322
    %n2s
    %nortsout
    speedleft=speed_y(w12_n);
    speedavg=speed_y((w12_n+w15_s)/2);
    speedrite=speed_y(w15_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F12_n+F15_s-speed*(w15_s-w12_n))*wgtnuv2/2;
    F_dql=(DF12_n+speed)*wgtnuv2/2;
    F_dqr=(DF15_s-speed)*wgtnuv2/2;
    residual_cell_312=residual_cell_312+psinT*Flux;
    residual_cell_322=residual_cell_322-psisT*Flux;
    Jac_cell_312=Jac_cell_312+psinT*F_dql*psin;
    Jac_nort_othr_312=Jac_nort_othr_312+psinT*F_dqr*psis;
    Jac_cell_322=Jac_cell_322-psisT*F_dqr*psis;
    Jac_sout_othr_322=Jac_sout_othr_322-psisT*F_dql*psin;

    %Interface15-18
    %cellleft322
    %cellrite332
    %n2s
    %nortsout
    speedleft=speed_y(w15_n);
    speedavg=speed_y((w15_n+w18_s)/2);
    speedrite=speed_y(w18_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F15_n+F18_s-speed*(w18_s-w15_n))*wgtnuv2/2;
    F_dql=(DF15_n+speed)*wgtnuv2/2;
    F_dqr=(DF18_s-speed)*wgtnuv2/2;
    residual_cell_322=residual_cell_322+psinT*Flux;
    residual_cell_332=residual_cell_332-psisT*Flux;
    Jac_cell_322=Jac_cell_322+psinT*F_dql*psin;
    Jac_nort_othr_322=Jac_nort_othr_322+psinT*F_dqr*psis;
    Jac_cell_332=Jac_cell_332-psisT*F_dqr*psis;
    Jac_sout_othr_332=Jac_sout_othr_332-psisT*F_dql*psin;

    %Interface19-22
    %cellleft113
    %cellrite123
    %n2s
    %nortsout
    speedleft=speed_y(w19_n);
    speedavg=speed_y((w19_n+w22_s)/2);
    speedrite=speed_y(w22_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F19_n+F22_s-speed*(w22_s-w19_n))*wgtnuv2/2;
    F_dql=(DF19_n+speed)*wgtnuv2/2;
    F_dqr=(DF22_s-speed)*wgtnuv2/2;
    residual_cell_113=residual_cell_113+psinT*Flux;
    residual_cell_123=residual_cell_123-psisT*Flux;
    Jac_cell_113=Jac_cell_113+psinT*F_dql*psin;
    Jac_nort_othr_113=Jac_nort_othr_113+psinT*F_dqr*psis;
    Jac_cell_123=Jac_cell_123-psisT*F_dqr*psis;
    Jac_sout_othr_123=Jac_sout_othr_123-psisT*F_dql*psin;

    %Interface22-25
    %cellleft123
    %cellrite133
    %n2s
    %nortsout
    speedleft=speed_y(w22_n);
    speedavg=speed_y((w22_n+w25_s)/2);
    speedrite=speed_y(w25_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F22_n+F25_s-speed*(w25_s-w22_n))*wgtnuv2/2;
    F_dql=(DF22_n+speed)*wgtnuv2/2;
    F_dqr=(DF25_s-speed)*wgtnuv2/2;
    residual_cell_123=residual_cell_123+psinT*Flux;
    residual_cell_133=residual_cell_133-psisT*Flux;
    Jac_cell_123=Jac_cell_123+psinT*F_dql*psin;
    Jac_nort_othr_123=Jac_nort_othr_123+psinT*F_dqr*psis;
    Jac_cell_133=Jac_cell_133-psisT*F_dqr*psis;
    Jac_sout_othr_133=Jac_sout_othr_133-psisT*F_dql*psin;

    %Interface20-23
    %cellleft213
    %cellrite223
    %n2s
    %nortsout
    speedleft=speed_y(w20_n);
    speedavg=speed_y((w20_n+w23_s)/2);
    speedrite=speed_y(w23_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F20_n+F23_s-speed*(w23_s-w20_n))*wgtnuv2/2;
    F_dql=(DF20_n+speed)*wgtnuv2/2;
    F_dqr=(DF23_s-speed)*wgtnuv2/2;
    residual_cell_213=residual_cell_213+psinT*Flux;
    residual_cell_223=residual_cell_223-psisT*Flux;
    Jac_cell_213=Jac_cell_213+psinT*F_dql*psin;
    Jac_nort_othr_213=Jac_nort_othr_213+psinT*F_dqr*psis;
    Jac_cell_223=Jac_cell_223-psisT*F_dqr*psis;
    Jac_sout_othr_223=Jac_sout_othr_223-psisT*F_dql*psin;

    %Interface23-26
    %cellleft223
    %cellrite233
    %n2s
    %nortsout
    speedleft=speed_y(w23_n);
    speedavg=speed_y((w23_n+w26_s)/2);
    speedrite=speed_y(w26_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F23_n+F26_s-speed*(w26_s-w23_n))*wgtnuv2/2;
    F_dql=(DF23_n+speed)*wgtnuv2/2;
    F_dqr=(DF26_s-speed)*wgtnuv2/2;
    residual_cell_223=residual_cell_223+psinT*Flux;
    residual_cell_233=residual_cell_233-psisT*Flux;
    Jac_cell_223=Jac_cell_223+psinT*F_dql*psin;
    Jac_nort_othr_223=Jac_nort_othr_223+psinT*F_dqr*psis;
    Jac_cell_233=Jac_cell_233-psisT*F_dqr*psis;
    Jac_sout_othr_233=Jac_sout_othr_233-psisT*F_dql*psin;

    %Interface21-24
    %cellleft313
    %cellrite323
    %n2s
    %nortsout
    speedleft=speed_y(w21_n);
    speedavg=speed_y((w21_n+w24_s)/2);
    speedrite=speed_y(w24_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F21_n+F24_s-speed*(w24_s-w21_n))*wgtnuv2/2;
    F_dql=(DF21_n+speed)*wgtnuv2/2;
    F_dqr=(DF24_s-speed)*wgtnuv2/2;
    residual_cell_313=residual_cell_313+psinT*Flux;
    residual_cell_323=residual_cell_323-psisT*Flux;
    Jac_cell_313=Jac_cell_313+psinT*F_dql*psin;
    Jac_nort_othr_313=Jac_nort_othr_313+psinT*F_dqr*psis;
    Jac_cell_323=Jac_cell_323-psisT*F_dqr*psis;
    Jac_sout_othr_323=Jac_sout_othr_323-psisT*F_dql*psin;

    %Interface24-27
    %cellleft323
    %cellrite333
    %n2s
    %nortsout
    speedleft=speed_y(w24_n);
    speedavg=speed_y((w24_n+w27_s)/2);
    speedrite=speed_y(w27_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F24_n+F27_s-speed*(w27_s-w24_n))*wgtnuv2/2;
    F_dql=(DF24_n+speed)*wgtnuv2/2;
    F_dqr=(DF27_s-speed)*wgtnuv2/2;
    residual_cell_323=residual_cell_323+psinT*Flux;
    residual_cell_333=residual_cell_333-psisT*Flux;
    Jac_cell_323=Jac_cell_323+psinT*F_dql*psin;
    Jac_nort_othr_323=Jac_nort_othr_323+psinT*F_dqr*psis;
    Jac_cell_333=Jac_cell_333-psisT*F_dqr*psis;
    Jac_sout_othr_333=Jac_sout_othr_333-psisT*F_dql*psin;


    %Interface1-10
    %cellleft111
    %cellrite112
    %u3d
    %upprdown
    speedleft=speed_z(w1_u);
    speedavg=speed_z((w1_u+w10_d)/2);
    speedrite=speed_z(w10_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F1_u+F10_d-speed*(w10_d-w1_u))*wgtnuv3/2;
    F_dql=(DF1_u+speed)*wgtnuv3/2;
    F_dqr=(DF10_d-speed)*wgtnuv3/2;
    residual_cell_111=residual_cell_111+psiuT*Flux;
    residual_cell_112=residual_cell_112-psidT*Flux;
    Jac_cell_111=Jac_cell_111+psiuT*F_dql*psiu;
    Jac_uppr_othr_111=Jac_uppr_othr_111+psiuT*F_dqr*psid;
    Jac_cell_112=Jac_cell_112-psidT*F_dqr*psid;
    Jac_down_othr_112=Jac_down_othr_112-psidT*F_dql*psiu;

    %Interface2-11
    %cellleft211
    %cellrite212
    %u3d
    %upprdown
    speedleft=speed_z(w2_u);
    speedavg=speed_z((w2_u+w11_d)/2);
    speedrite=speed_z(w11_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F2_u+F11_d-speed*(w11_d-w2_u))*wgtnuv3/2;
    F_dql=(DF2_u+speed)*wgtnuv3/2;
    F_dqr=(DF11_d-speed)*wgtnuv3/2;
    residual_cell_211=residual_cell_211+psiuT*Flux;
    residual_cell_212=residual_cell_212-psidT*Flux;
    Jac_cell_211=Jac_cell_211+psiuT*F_dql*psiu;
    Jac_uppr_othr_211=Jac_uppr_othr_211+psiuT*F_dqr*psid;
    Jac_cell_212=Jac_cell_212-psidT*F_dqr*psid;
    Jac_down_othr_212=Jac_down_othr_212-psidT*F_dql*psiu;

    %Interface3-12
    %cellleft311
    %cellrite312
    %u3d
    %upprdown
    speedleft=speed_z(w3_u);
    speedavg=speed_z((w3_u+w12_d)/2);
    speedrite=speed_z(w12_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F3_u+F12_d-speed*(w12_d-w3_u))*wgtnuv3/2;
    F_dql=(DF3_u+speed)*wgtnuv3/2;
    F_dqr=(DF12_d-speed)*wgtnuv3/2;
    residual_cell_311=residual_cell_311+psiuT*Flux;
    residual_cell_312=residual_cell_312-psidT*Flux;
    Jac_cell_311=Jac_cell_311+psiuT*F_dql*psiu;
    Jac_uppr_othr_311=Jac_uppr_othr_311+psiuT*F_dqr*psid;
    Jac_cell_312=Jac_cell_312-psidT*F_dqr*psid;
    Jac_down_othr_312=Jac_down_othr_312-psidT*F_dql*psiu;

    %Interface4-13
    %cellleft121
    %cellrite122
    %u3d
    %upprdown
    speedleft=speed_z(w4_u);
    speedavg=speed_z((w4_u+w13_d)/2);
    speedrite=speed_z(w13_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F4_u+F13_d-speed*(w13_d-w4_u))*wgtnuv3/2;
    F_dql=(DF4_u+speed)*wgtnuv3/2;
    F_dqr=(DF13_d-speed)*wgtnuv3/2;
    residual_cell_121=residual_cell_121+psiuT*Flux;
    residual_cell_122=residual_cell_122-psidT*Flux;
    Jac_cell_121=Jac_cell_121+psiuT*F_dql*psiu;
    Jac_uppr_othr_121=Jac_uppr_othr_121+psiuT*F_dqr*psid;
    Jac_cell_122=Jac_cell_122-psidT*F_dqr*psid;
    Jac_down_othr_122=Jac_down_othr_122-psidT*F_dql*psiu;

    %Interface5-14
    %cellleft221
    %cellrite222
    %u3d
    %upprdown
    speedleft=speed_z(w5_u);
    speedavg=speed_z((w5_u+w14_d)/2);
    speedrite=speed_z(w14_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F5_u+F14_d-speed*(w14_d-w5_u))*wgtnuv3/2;
    F_dql=(DF5_u+speed)*wgtnuv3/2;
    F_dqr=(DF14_d-speed)*wgtnuv3/2;
    residual_cell_221=residual_cell_221+psiuT*Flux;
    residual_cell_222=residual_cell_222-psidT*Flux;
    Jac_cell_221=Jac_cell_221+psiuT*F_dql*psiu;
    Jac_uppr_othr_221=Jac_uppr_othr_221+psiuT*F_dqr*psid;
    Jac_cell_222=Jac_cell_222-psidT*F_dqr*psid;
    Jac_down_othr_222=Jac_down_othr_222-psidT*F_dql*psiu;

    %Interface6-15
    %cellleft321
    %cellrite322
    %u3d
    %upprdown
    speedleft=speed_z(w6_u);
    speedavg=speed_z((w6_u+w15_d)/2);
    speedrite=speed_z(w15_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F6_u+F15_d-speed*(w15_d-w6_u))*wgtnuv3/2;
    F_dql=(DF6_u+speed)*wgtnuv3/2;
    F_dqr=(DF15_d-speed)*wgtnuv3/2;
    residual_cell_321=residual_cell_321+psiuT*Flux;
    residual_cell_322=residual_cell_322-psidT*Flux;
    Jac_cell_321=Jac_cell_321+psiuT*F_dql*psiu;
    Jac_uppr_othr_321=Jac_uppr_othr_321+psiuT*F_dqr*psid;
    Jac_cell_322=Jac_cell_322-psidT*F_dqr*psid;
    Jac_down_othr_322=Jac_down_othr_322-psidT*F_dql*psiu;

    %Interface7-16
    %cellleft131
    %cellrite132
    %u3d
    %upprdown
    speedleft=speed_z(w7_u);
    speedavg=speed_z((w7_u+w16_d)/2);
    speedrite=speed_z(w16_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F7_u+F16_d-speed*(w16_d-w7_u))*wgtnuv3/2;
    F_dql=(DF7_u+speed)*wgtnuv3/2;
    F_dqr=(DF16_d-speed)*wgtnuv3/2;
    residual_cell_131=residual_cell_131+psiuT*Flux;
    residual_cell_132=residual_cell_132-psidT*Flux;
    Jac_cell_131=Jac_cell_131+psiuT*F_dql*psiu;
    Jac_uppr_othr_131=Jac_uppr_othr_131+psiuT*F_dqr*psid;
    Jac_cell_132=Jac_cell_132-psidT*F_dqr*psid;
    Jac_down_othr_132=Jac_down_othr_132-psidT*F_dql*psiu;

    %Interface8-17
    %cellleft231
    %cellrite232
    %u3d
    %upprdown
    speedleft=speed_z(w8_u);
    speedavg=speed_z((w8_u+w17_d)/2);
    speedrite=speed_z(w17_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F8_u+F17_d-speed*(w17_d-w8_u))*wgtnuv3/2;
    F_dql=(DF8_u+speed)*wgtnuv3/2;
    F_dqr=(DF17_d-speed)*wgtnuv3/2;
    residual_cell_231=residual_cell_231+psiuT*Flux;
    residual_cell_232=residual_cell_232-psidT*Flux;
    Jac_cell_231=Jac_cell_231+psiuT*F_dql*psiu;
    Jac_uppr_othr_231=Jac_uppr_othr_231+psiuT*F_dqr*psid;
    Jac_cell_232=Jac_cell_232-psidT*F_dqr*psid;
    Jac_down_othr_232=Jac_down_othr_232-psidT*F_dql*psiu;

    %Interface9-18
    %cellleft331
    %cellrite332
    %u3d
    %upprdown
    speedleft=speed_z(w9_u);
    speedavg=speed_z((w9_u+w18_d)/2);
    speedrite=speed_z(w18_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F9_u+F18_d-speed*(w18_d-w9_u))*wgtnuv3/2;
    F_dql=(DF9_u+speed)*wgtnuv3/2;
    F_dqr=(DF18_d-speed)*wgtnuv3/2;
    residual_cell_331=residual_cell_331+psiuT*Flux;
    residual_cell_332=residual_cell_332-psidT*Flux;
    Jac_cell_331=Jac_cell_331+psiuT*F_dql*psiu;
    Jac_uppr_othr_331=Jac_uppr_othr_331+psiuT*F_dqr*psid;
    Jac_cell_332=Jac_cell_332-psidT*F_dqr*psid;
    Jac_down_othr_332=Jac_down_othr_332-psidT*F_dql*psiu;

    %Interface10-19
    %cellleft112
    %cellrite113
    %u3d
    %upprdown
    speedleft=speed_z(w10_u);
    speedavg=speed_z((w10_u+w19_d)/2);
    speedrite=speed_z(w19_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F10_u+F19_d-speed*(w19_d-w10_u))*wgtnuv3/2;
    F_dql=(DF10_u+speed)*wgtnuv3/2;
    F_dqr=(DF19_d-speed)*wgtnuv3/2;
    residual_cell_112=residual_cell_112+psiuT*Flux;
    residual_cell_113=residual_cell_113-psidT*Flux;
    Jac_cell_112=Jac_cell_112+psiuT*F_dql*psiu;
    Jac_uppr_othr_112=Jac_uppr_othr_112+psiuT*F_dqr*psid;
    Jac_cell_113=Jac_cell_113-psidT*F_dqr*psid;
    Jac_down_othr_113=Jac_down_othr_113-psidT*F_dql*psiu;

    %Interface11-20
    %cellleft212
    %cellrite213
    %u3d
    %upprdown
    speedleft=speed_z(w11_u);
    speedavg=speed_z((w11_u+w20_d)/2);
    speedrite=speed_z(w20_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F11_u+F20_d-speed*(w20_d-w11_u))*wgtnuv3/2;
    F_dql=(DF11_u+speed)*wgtnuv3/2;
    F_dqr=(DF20_d-speed)*wgtnuv3/2;
    residual_cell_212=residual_cell_212+psiuT*Flux;
    residual_cell_213=residual_cell_213-psidT*Flux;
    Jac_cell_212=Jac_cell_212+psiuT*F_dql*psiu;
    Jac_uppr_othr_212=Jac_uppr_othr_212+psiuT*F_dqr*psid;
    Jac_cell_213=Jac_cell_213-psidT*F_dqr*psid;
    Jac_down_othr_213=Jac_down_othr_213-psidT*F_dql*psiu;

    %Interface12-21
    %cellleft312
    %cellrite313
    %u3d
    %upprdown
    speedleft=speed_z(w12_u);
    speedavg=speed_z((w12_u+w21_d)/2);
    speedrite=speed_z(w21_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F12_u+F21_d-speed*(w21_d-w12_u))*wgtnuv3/2;
    F_dql=(DF12_u+speed)*wgtnuv3/2;
    F_dqr=(DF21_d-speed)*wgtnuv3/2;
    residual_cell_312=residual_cell_312+psiuT*Flux;
    residual_cell_313=residual_cell_313-psidT*Flux;
    Jac_cell_312=Jac_cell_312+psiuT*F_dql*psiu;
    Jac_uppr_othr_312=Jac_uppr_othr_312+psiuT*F_dqr*psid;
    Jac_cell_313=Jac_cell_313-psidT*F_dqr*psid;
    Jac_down_othr_313=Jac_down_othr_313-psidT*F_dql*psiu;

    %Interface13-22
    %cellleft122
    %cellrite123
    %u3d
    %upprdown
    speedleft=speed_z(w13_u);
    speedavg=speed_z((w13_u+w22_d)/2);
    speedrite=speed_z(w22_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F13_u+F22_d-speed*(w22_d-w13_u))*wgtnuv3/2;
    F_dql=(DF13_u+speed)*wgtnuv3/2;
    F_dqr=(DF22_d-speed)*wgtnuv3/2;
    residual_cell_122=residual_cell_122+psiuT*Flux;
    residual_cell_123=residual_cell_123-psidT*Flux;
    Jac_cell_122=Jac_cell_122+psiuT*F_dql*psiu;
    Jac_uppr_othr_122=Jac_uppr_othr_122+psiuT*F_dqr*psid;
    Jac_cell_123=Jac_cell_123-psidT*F_dqr*psid;
    Jac_down_othr_123=Jac_down_othr_123-psidT*F_dql*psiu;

    %Interface14-23
    %cellleft222
    %cellrite223
    %u3d
    %upprdown
    speedleft=speed_z(w14_u);
    speedavg=speed_z((w14_u+w23_d)/2);
    speedrite=speed_z(w23_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F14_u+F23_d-speed*(w23_d-w14_u))*wgtnuv3/2;
    F_dql=(DF14_u+speed)*wgtnuv3/2;
    F_dqr=(DF23_d-speed)*wgtnuv3/2;
    residual_cell_222=residual_cell_222+psiuT*Flux;
    residual_cell_223=residual_cell_223-psidT*Flux;
    Jac_cell_222=Jac_cell_222+psiuT*F_dql*psiu;
    Jac_uppr_othr_222=Jac_uppr_othr_222+psiuT*F_dqr*psid;
    Jac_cell_223=Jac_cell_223-psidT*F_dqr*psid;
    Jac_down_othr_223=Jac_down_othr_223-psidT*F_dql*psiu;

    %Interface15-24
    %cellleft322
    %cellrite323
    %u3d
    %upprdown
    speedleft=speed_z(w15_u);
    speedavg=speed_z((w15_u+w24_d)/2);
    speedrite=speed_z(w24_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F15_u+F24_d-speed*(w24_d-w15_u))*wgtnuv3/2;
    F_dql=(DF15_u+speed)*wgtnuv3/2;
    F_dqr=(DF24_d-speed)*wgtnuv3/2;
    residual_cell_322=residual_cell_322+psiuT*Flux;
    residual_cell_323=residual_cell_323-psidT*Flux;
    Jac_cell_322=Jac_cell_322+psiuT*F_dql*psiu;
    Jac_uppr_othr_322=Jac_uppr_othr_322+psiuT*F_dqr*psid;
    Jac_cell_323=Jac_cell_323-psidT*F_dqr*psid;
    Jac_down_othr_323=Jac_down_othr_323-psidT*F_dql*psiu;

    %Interface16-25
    %cellleft132
    %cellrite133
    %u3d
    %upprdown
    speedleft=speed_z(w16_u);
    speedavg=speed_z((w16_u+w25_d)/2);
    speedrite=speed_z(w25_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F16_u+F25_d-speed*(w25_d-w16_u))*wgtnuv3/2;
    F_dql=(DF16_u+speed)*wgtnuv3/2;
    F_dqr=(DF25_d-speed)*wgtnuv3/2;
    residual_cell_132=residual_cell_132+psiuT*Flux;
    residual_cell_133=residual_cell_133-psidT*Flux;
    Jac_cell_132=Jac_cell_132+psiuT*F_dql*psiu;
    Jac_uppr_othr_132=Jac_uppr_othr_132+psiuT*F_dqr*psid;
    Jac_cell_133=Jac_cell_133-psidT*F_dqr*psid;
    Jac_down_othr_133=Jac_down_othr_133-psidT*F_dql*psiu;

    %Interface17-26
    %cellleft232
    %cellrite233
    %u3d
    %upprdown
    speedleft=speed_z(w17_u);
    speedavg=speed_z((w17_u+w26_d)/2);
    speedrite=speed_z(w26_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F17_u+F26_d-speed*(w26_d-w17_u))*wgtnuv3/2;
    F_dql=(DF17_u+speed)*wgtnuv3/2;
    F_dqr=(DF26_d-speed)*wgtnuv3/2;
    residual_cell_232=residual_cell_232+psiuT*Flux;
    residual_cell_233=residual_cell_233-psidT*Flux;
    Jac_cell_232=Jac_cell_232+psiuT*F_dql*psiu;
    Jac_uppr_othr_232=Jac_uppr_othr_232+psiuT*F_dqr*psid;
    Jac_cell_233=Jac_cell_233-psidT*F_dqr*psid;
    Jac_down_othr_233=Jac_down_othr_233-psidT*F_dql*psiu;

    %Interface18-27
    %cellleft332
    %cellrite333
    %u3d
    %upprdown
    speedleft=speed_z(w18_u);
    speedavg=speed_z((w18_u+w27_d)/2);
    speedrite=speed_z(w27_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F18_u+F27_d-speed*(w27_d-w18_u))*wgtnuv3/2;
    F_dql=(DF18_u+speed)*wgtnuv3/2;
    F_dqr=(DF27_d-speed)*wgtnuv3/2;
    residual_cell_332=residual_cell_332+psiuT*Flux;
    residual_cell_333=residual_cell_333-psidT*Flux;
    Jac_cell_332=Jac_cell_332+psiuT*F_dql*psiu;
    Jac_uppr_othr_332=Jac_uppr_othr_332+psiuT*F_dqr*psid;
    Jac_cell_333=Jac_cell_333-psidT*F_dqr*psid;
    Jac_down_othr_333=Jac_down_othr_333-psidT*F_dql*psiu;
    
    
    
   
    
end

indices = @(i) (1:thetaT) + (i-1)*thetaT;
                               
residual(indices(1),1)=residual_cell_111;%
residual(indices(2),1)=residual_cell_211;%
residual(indices(3),1)=residual_cell_311;%
residual(indices(4),1)=residual_cell_121;%
residual(indices(5),1)=residual_cell_221;%
residual(indices(6),1)=residual_cell_321;%
residual(indices(7),1)=residual_cell_131;%
residual(indices(8),1)=residual_cell_231;%
residual(indices(9),1)=residual_cell_331;%
residual(indices(10),1)=residual_cell_112;%
residual(indices(11),1)=residual_cell_212;%
residual(indices(12),1)=residual_cell_312;%
residual(indices(13),1)=residual_cell_122;%
residual(indices(14),1)=residual_cell_222;%
residual(indices(15),1)=residual_cell_322;%
residual(indices(16),1)=residual_cell_132;%
residual(indices(17),1)=residual_cell_232;%
residual(indices(18),1)=residual_cell_332;%
residual(indices(19),1)=residual_cell_113;%
residual(indices(20),1)=residual_cell_213;%
residual(indices(21),1)=residual_cell_313;%
residual(indices(22),1)=residual_cell_123;%
residual(indices(23),1)=residual_cell_223;%
residual(indices(24),1)=residual_cell_323;%
residual(indices(25),1)=residual_cell_133;%
residual(indices(26),1)=residual_cell_233;%
residual(indices(27),1)=residual_cell_333;%

Jacobian(indices(1),incides(1))=Jacobian_cell_111+Jac_east_flux_111+Jac_west_trunc_111+Jac_nort_flux_111+Jac_sout_trunc_111+Jac_uppr_flux_111+Jac_down_trunc_111;%
Jacobian(indices(2),incides(2))=Jacobian_cell_211+Jac_east_flux_211+Jac_west_flux_211+Jac_nort_flux_211+Jac_sout_trunc_211+Jac_uppr_flux_211+Jac_down_trunc_211;%
Jacobian(indices(3),incides(3))=Jacobian_cell_311+Jac_east_trunc_311+Jac_west_flux_311+Jac_nort_flux_311+Jac_sout_trunc_311+Jac_uppr_flux_311+Jac_down_trunc_311;%
Jacobian(indices(4),incides(4))=Jacobian_cell_121+Jac_east_flux_121+Jac_west_trunc_121+Jac_nort_flux_121+Jac_sout_flux_121+Jac_uppr_flux_121+Jac_down_trunc_121;%
Jacobian(indices(5),incides(5))=Jacobian_cell_221+Jac_east_flux_221+Jac_west_flux_221+Jac_nort_flux_221+Jac_sout_flux_221+Jac_uppr_flux_221+Jac_down_trunc_221;%
Jacobian(indices(6),incides(6))=Jacobian_cell_321+Jac_east_trunc_321+Jac_west_flux_321+Jac_nort_flux_321+Jac_sout_flux_321+Jac_uppr_flux_321+Jac_down_trunc_321;%
Jacobian(indices(7),incides(7))=Jacobian_cell_131+Jac_east_flux_131+Jac_west_trunc_131+Jac_nort_trunc_131+Jac_sout_flux_131+Jac_uppr_flux_131+Jac_down_trunc_131;%
Jacobian(indices(8),incides(8))=Jacobian_cell_231+Jac_east_flux_231+Jac_west_flux_231+Jac_nort_trunc_231+Jac_sout_flux_231+Jac_uppr_flux_231+Jac_down_trunc_231;%
Jacobian(indices(9),incides(9))=Jacobian_cell_331+Jac_east_trunc_331+Jac_west_flux_331+Jac_nort_trunc_331+Jac_sout_flux_331+Jac_uppr_flux_331+Jac_down_trunc_331;%
Jacobian(indices(10),incides(10))=Jacobian_cell_112+Jac_east_flux_112+Jac_west_trunc_112+Jac_nort_flux_112+Jac_sout_trunc_112+Jac_uppr_flux_112+Jac_down_flux_112;%
Jacobian(indices(11),incides(11))=Jacobian_cell_212+Jac_east_flux_212+Jac_west_flux_212+Jac_nort_flux_212+Jac_sout_trunc_212+Jac_uppr_flux_212+Jac_down_flux_212;%
Jacobian(indices(12),incides(12))=Jacobian_cell_312+Jac_east_trunc_312+Jac_west_flux_312+Jac_nort_flux_312+Jac_sout_trunc_312+Jac_uppr_flux_312+Jac_down_flux_312;%
Jacobian(indices(13),incides(13))=Jacobian_cell_122+Jac_east_flux_122+Jac_west_trunc_122+Jac_nort_flux_122+Jac_sout_flux_122+Jac_uppr_flux_122+Jac_down_flux_122;%
Jacobian(indices(14),incides(14))=Jacobian_cell_222+Jac_east_flux_222+Jac_west_flux_222+Jac_nort_flux_222+Jac_sout_flux_222+Jac_uppr_flux_222+Jac_down_flux_222;%
Jacobian(indices(15),incides(15))=Jacobian_cell_322+Jac_east_trunc_322+Jac_west_flux_322+Jac_nort_flux_322+Jac_sout_flux_322+Jac_uppr_flux_322+Jac_down_flux_322;%
Jacobian(indices(16),incides(16))=Jacobian_cell_132+Jac_east_flux_132+Jac_west_trunc_132+Jac_nort_trunc_132+Jac_sout_flux_132+Jac_uppr_flux_132+Jac_down_flux_132;%
Jacobian(indices(17),incides(17))=Jacobian_cell_232+Jac_east_flux_232+Jac_west_flux_232+Jac_nort_trunc_232+Jac_sout_flux_232+Jac_uppr_flux_232+Jac_down_flux_232;%
Jacobian(indices(18),incides(18))=Jacobian_cell_332+Jac_east_trunc_332+Jac_west_flux_332+Jac_nort_trunc_332+Jac_sout_flux_332+Jac_uppr_flux_332+Jac_down_flux_332;%
Jacobian(indices(19),incides(19))=Jacobian_cell_113+Jac_east_flux_113+Jac_west_trunc_113+Jac_nort_flux_113+Jac_sout_trunc_113+Jac_uppr_trunc_113+Jac_down_flux_113;%
Jacobian(indices(20),incides(20))=Jacobian_cell_213+Jac_east_flux_213+Jac_west_flux_213+Jac_nort_flux_213+Jac_sout_trunc_213+Jac_uppr_trunc_213+Jac_down_flux_213;%
Jacobian(indices(21),incides(21))=Jacobian_cell_313+Jac_east_trunc_313+Jac_west_flux_313+Jac_nort_flux_313+Jac_sout_trunc_313+Jac_uppr_trunc_313+Jac_down_flux_313;%
Jacobian(indices(22),incides(22))=Jacobian_cell_123+Jac_east_flux_123+Jac_west_trunc_123+Jac_nort_flux_123+Jac_sout_flux_123+Jac_uppr_trunc_123+Jac_down_flux_123;%
Jacobian(indices(23),incides(23))=Jacobian_cell_223+Jac_east_flux_223+Jac_west_flux_223+Jac_nort_flux_223+Jac_sout_flux_223+Jac_uppr_trunc_223+Jac_down_flux_223;%
Jacobian(indices(24),incides(24))=Jacobian_cell_323+Jac_east_trunc_323+Jac_west_flux_323+Jac_nort_flux_323+Jac_sout_flux_323+Jac_uppr_trunc_323+Jac_down_flux_323;%
Jacobian(indices(25),incides(25))=Jacobian_cell_133+Jac_east_flux_133+Jac_west_trunc_133+Jac_nort_trunc_133+Jac_sout_flux_133+Jac_uppr_trunc_133+Jac_down_flux_133;%
Jacobian(indices(26),incides(26))=Jacobian_cell_233+Jac_east_flux_233+Jac_west_flux_233+Jac_nort_trunc_233+Jac_sout_flux_233+Jac_uppr_trunc_233+Jac_down_flux_233;%
Jacobian(indices(27),incides(27))=Jacobian_cell_333+Jac_east_trunc_333+Jac_west_flux_333+Jac_nort_trunc_333+Jac_sout_flux_333+Jac_uppr_trunc_333+Jac_down_flux_333;%

Jacobian(indices(1),incides(2)) = Jac_east_othr_111;%
Jacobian(indices(2),incides(3)) = Jac_east_othr_211;%
%Jacobian(indices(3),incides(4)) = Jac_east_othr_311;%
Jacobian(indices(4),incides(5)) = Jac_east_othr_121;%
Jacobian(indices(5),incides(6)) = Jac_east_othr_221;%
%Jacobian(indices(6),incides(7)) = Jac_east_othr_321;%
Jacobian(indices(7),incides(8)) = Jac_east_othr_131;%
Jacobian(indices(8),incides(9)) = Jac_east_othr_231;%
%Jacobian(indices(9),incides(10)) = Jac_east_othr_331;%
Jacobian(indices(10),incides(11)) = Jac_east_othr_112;%
Jacobian(indices(11),incides(12)) = Jac_east_othr_212;%
%Jacobian(indices(12),incides(13)) = Jac_east_othr_312;%
Jacobian(indices(13),incides(14)) = Jac_east_othr_122;%
Jacobian(indices(14),incides(15)) = Jac_east_othr_222;%
%Jacobian(indices(15),incides(16)) = Jac_east_othr_322;%
Jacobian(indices(16),incides(17)) = Jac_east_othr_132;%
Jacobian(indices(17),incides(18)) = Jac_east_othr_232;%
%Jacobian(indices(18),incides(19)) = Jac_east_othr_332;%
Jacobian(indices(19),incides(20)) = Jac_east_othr_113;%
Jacobian(indices(20),incides(21)) = Jac_east_othr_213;%
%Jacobian(indices(21),incides(22)) = Jac_east_othr_313;%
Jacobian(indices(22),incides(23)) = Jac_east_othr_123;%
Jacobian(indices(23),incides(24)) = Jac_east_othr_223;%
%Jacobian(indices(24),incides(25)) = Jac_east_othr_323;%
Jacobian(indices(25),incides(26)) = Jac_east_othr_133;%
Jacobian(indices(26),incides(27)) = Jac_east_othr_233;%
%;%
%Jacobian(indices(1),incides(0)) = Jac_west_othr_111;%
Jacobian(indices(2),incides(1)) = Jac_west_othr_211;%
Jacobian(indices(3),incides(2)) = Jac_west_othr_311;%
%Jacobian(indices(4),incides(3)) = Jac_west_othr_121;%
Jacobian(indices(5),incides(4)) = Jac_west_othr_221;%
Jacobian(indices(6),incides(5)) = Jac_west_othr_321;%
%Jacobian(indices(7),incides(6)) = Jac_west_othr_131;%
Jacobian(indices(8),incides(7)) = Jac_west_othr_231;%
Jacobian(indices(9),incides(8)) = Jac_west_othr_331;%
%Jacobian(indices(10),incides(9)) = Jac_west_othr_112;%
Jacobian(indices(11),incides(10)) = Jac_west_othr_212;%
Jacobian(indices(12),incides(11)) = Jac_west_othr_312;%
%Jacobian(indices(13),incides(12)) = Jac_west_othr_122;%
Jacobian(indices(14),incides(13)) = Jac_west_othr_222;%
Jacobian(indices(15),incides(14)) = Jac_west_othr_322;%
%Jacobian(indices(16),incides(15)) = Jac_west_othr_132;%
Jacobian(indices(17),incides(16)) = Jac_west_othr_232;%
Jacobian(indices(18),incides(17)) = Jac_west_othr_332;%
%Jacobian(indices(19),incides(18)) = Jac_west_othr_113;%
Jacobian(indices(20),incides(19)) = Jac_west_othr_213;%
Jacobian(indices(21),incides(20)) = Jac_west_othr_313;%
%Jacobian(indices(22),incides(21)) = Jac_west_othr_123;%
Jacobian(indices(23),incides(22)) = Jac_west_othr_223;%
Jacobian(indices(24),incides(23)) = Jac_west_othr_323;%
%Jacobian(indices(25),incides(24)) = Jac_west_othr_133;%
Jacobian(indices(26),incides(25)) = Jac_west_othr_233;%
Jacobian(indices(27),incides(26)) = Jac_west_othr_333;%
%;%
Jacobian(indices(1),incides(4)) = Jac_nort_othr_111;%
Jacobian(indices(2),incides(5)) = Jac_nort_othr_211;%
Jacobian(indices(3),incides(6)) = Jac_nort_othr_311;%
Jacobian(indices(4),incides(7)) = Jac_nort_othr_121;%
Jacobian(indices(5),incides(8)) = Jac_nort_othr_221;%
Jacobian(indices(6),incides(9)) = Jac_nort_othr_321;%
%Jacobian(indices(7),incides(10)) = Jac_nort_othr_131;%
%Jacobian(indices(8),incides(11)) = Jac_nort_othr_231;%
%Jacobian(indices(9),incides(12)) = Jac_nort_othr_331;%
Jacobian(indices(10),incides(13)) = Jac_nort_othr_112;%
Jacobian(indices(11),incides(14)) = Jac_nort_othr_212;%
Jacobian(indices(12),incides(15)) = Jac_nort_othr_312;%
Jacobian(indices(13),incides(16)) = Jac_nort_othr_122;%
Jacobian(indices(14),incides(17)) = Jac_nort_othr_222;%
Jacobian(indices(15),incides(18)) = Jac_nort_othr_322;%
%Jacobian(indices(16),incides(19)) = Jac_nort_othr_132;%
%Jacobian(indices(17),incides(20)) = Jac_nort_othr_232;%
%Jacobian(indices(18),incides(21)) = Jac_nort_othr_332;%
Jacobian(indices(19),incides(22)) = Jac_nort_othr_113;%
Jacobian(indices(20),incides(23)) = Jac_nort_othr_213;%
Jacobian(indices(21),incides(24)) = Jac_nort_othr_313;%
Jacobian(indices(22),incides(25)) = Jac_nort_othr_123;%
Jacobian(indices(23),incides(26)) = Jac_nort_othr_223;%
Jacobian(indices(24),incides(27)) = Jac_nort_othr_323;%
%Jacobian(indices(25),incides(28)) = Jac_nort_othr_133;%
%Jacobian(indices(26),incides(29)) = Jac_nort_othr_233;%
%;%
%Jacobian(indices(1),incides(-15)) = Jac_sout_othr_111;%
%Jacobian(indices(2),incides(-1)) = Jac_sout_othr_211;%
%Jacobian(indices(3),incides(0)) = Jac_sout_othr_311;%
%Jacobian(indices(4),incides(1)) = Jac_sout_othr_121;%
Jacobian(indices(5),incides(2)) = Jac_sout_othr_221;%
Jacobian(indices(6),incides(3)) = Jac_sout_othr_321;%
Jacobian(indices(7),incides(4)) = Jac_sout_othr_131;%
Jacobian(indices(8),incides(5)) = Jac_sout_othr_231;%
Jacobian(indices(9),incides(6)) = Jac_sout_othr_331;%
Jacobian(indices(10),incides(7)) = Jac_sout_othr_112;%
%Jacobian(indices(11),incides(8)) = Jac_sout_othr_212;%
%Jacobian(indices(12),incides(9)) = Jac_sout_othr_312;%
%Jacobian(indices(13),incides(10)) = Jac_sout_othr_122;%
Jacobian(indices(14),incides(11)) = Jac_sout_othr_222;%
Jacobian(indices(15),incides(12)) = Jac_sout_othr_322;%
Jacobian(indices(16),incides(13)) = Jac_sout_othr_132;%
Jacobian(indices(17),incides(14)) = Jac_sout_othr_232;%
Jacobian(indices(18),incides(15)) = Jac_sout_othr_332;%
Jacobian(indices(19),incides(16)) = Jac_sout_othr_113;%
%Jacobian(indices(20),incides(17)) = Jac_sout_othr_213;%
%Jacobian(indices(21),incides(18)) = Jac_sout_othr_313;%
%Jacobian(indices(22),incides(19)) = Jac_sout_othr_123;%
Jacobian(indices(23),incides(20)) = Jac_sout_othr_223;%
Jacobian(indices(24),incides(21)) = Jac_sout_othr_323;%
Jacobian(indices(25),incides(22)) = Jac_sout_othr_133;%
Jacobian(indices(26),incides(23)) = Jac_sout_othr_233;%
Jacobian(indices(27),incides(24)) = Jac_sout_othr_333;%
%;%
Jacobian(indices(1),incides(10)) = Jac_uppr_othr_111;%
Jacobian(indices(2),incides(11)) = Jac_uppr_othr_211;%
Jacobian(indices(3),incides(12)) = Jac_uppr_othr_311;%
Jacobian(indices(4),incides(13)) = Jac_uppr_othr_121;%
Jacobian(indices(5),incides(14)) = Jac_uppr_othr_221;%
Jacobian(indices(6),incides(15)) = Jac_uppr_othr_321;%
Jacobian(indices(7),incides(16)) = Jac_uppr_othr_131;%
Jacobian(indices(8),incides(17)) = Jac_uppr_othr_231;%
Jacobian(indices(9),incides(18)) = Jac_uppr_othr_331;%
Jacobian(indices(10),incides(19)) = Jac_uppr_othr_112;%
Jacobian(indices(11),incides(20)) = Jac_uppr_othr_212;%
Jacobian(indices(12),incides(21)) = Jac_uppr_othr_312;%
Jacobian(indices(13),incides(22)) = Jac_uppr_othr_122;%
Jacobian(indices(14),incides(23)) = Jac_uppr_othr_222;%
Jacobian(indices(15),incides(24)) = Jac_uppr_othr_322;%
Jacobian(indices(16),incides(25)) = Jac_uppr_othr_132;%
Jacobian(indices(17),incides(26)) = Jac_uppr_othr_232;%
Jacobian(indices(18),incides(27)) = Jac_uppr_othr_332;%
%Jacobian(indices(19),incides(28)) = Jac_uppr_othr_113;%
%Jacobian(indices(20),incides(29)) = Jac_uppr_othr_213;%
%Jacobian(indices(21),incides(30)) = Jac_uppr_othr_313;%
%Jacobian(indices(22),incides(31)) = Jac_uppr_othr_123;%
%Jacobian(indices(23),incides(32)) = Jac_uppr_othr_223;%
%Jacobian(indices(24),incides(33)) = Jac_uppr_othr_323;%
%Jacobian(indices(25),incides(34)) = Jac_uppr_othr_133;%
%Jacobian(indices(26),incides(35)) = Jac_uppr_othr_233;%
%;%
%Jacobian(indices(1),incides(-8)) = Jac_down_othr_111;%
%Jacobian(indices(2),incides(-7)) = Jac_down_othr_211;%
%Jacobian(indices(3),incides(-6)) = Jac_down_othr_311;%
%Jacobian(indices(4),incides(-5)) = Jac_down_othr_121;%
%Jacobian(indices(5),incides(-4)) = Jac_down_othr_221;%
%Jacobian(indices(6),incides(-3)) = Jac_down_othr_321;%
%Jacobian(indices(7),incides(-2)) = Jac_down_othr_131;%
%Jacobian(indices(8),incides(-1)) = Jac_down_othr_231;%
%Jacobian(indices(9),incides(0)) = Jac_down_othr_331;%
Jacobian(indices(10),incides(1)) = Jac_down_othr_112;%
Jacobian(indices(11),incides(2)) = Jac_down_othr_212;%
Jacobian(indices(12),incides(3)) = Jac_down_othr_312;%
Jacobian(indices(13),incides(4)) = Jac_down_othr_122;%
Jacobian(indices(14),incides(5)) = Jac_down_othr_222;%
Jacobian(indices(15),incides(6)) = Jac_down_othr_322;%
Jacobian(indices(16),incides(7)) = Jac_down_othr_132;%
Jacobian(indices(17),incides(8)) = Jac_down_othr_232;%
Jacobian(indices(18),incides(9)) = Jac_down_othr_332;%
Jacobian(indices(19),incides(10)) = Jac_down_othr_113;%
Jacobian(indices(20),incides(11)) = Jac_down_othr_213;%
Jacobian(indices(21),incides(12)) = Jac_down_othr_313;%
Jacobian(indices(22),incides(13)) = Jac_down_othr_123;%
Jacobian(indices(23),incides(14)) = Jac_down_othr_223;%
Jacobian(indices(24),incides(15)) = Jac_down_othr_323;%
Jacobian(indices(25),incides(16)) = Jac_down_othr_133;%
Jacobian(indices(26),incides(17)) = Jac_down_othr_233;%
Jacobian(indices(27),incides(18)) = Jac_down_othr_333;%
                                



end



