function [Jacobian,residual,auxiliary] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data,cellindex_center)
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
Dquadwgts = data.Dquadwgts;

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

vectphi_west = data.vectphi_west;
vectphi_east = data.vectphi_east;
vectphi_sout = data.vectphi_sout;
vectphi_nort = data.vectphi_nort;
vectphi_down = data.vectphi_down;
vectphi_uppr = data.vectphi_uppr;


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

qpast01 = DGregion_past(:,rx-1,ry-1,rz-1);
qpast02 = DGregion_past(:,rx+0,ry-1,rz-1);
qpast03 = DGregion_past(:,rx+1,ry-1,rz-1);

qpast04 = DGregion_past(:,rx-1,ry+0,rz-1);
qpast05 = DGregion_past(:,rx+0,ry+0,rz-1);
qpast06 = DGregion_past(:,rx+1,ry+0,rz-1);

qpast07 = DGregion_past(:,rx-1,ry+1,rz-1);
qpast08 = DGregion_past(:,rx+0,ry+1,rz-1);
qpast09 = DGregion_past(:,rx+1,ry+1,rz-1);

qpast10 = DGregion_past(:,rx-1,ry-1,rz);
qpast11 = DGregion_past(:,rx+0,ry-1,rz);
qpast12 = DGregion_past(:,rx+1,ry-1,rz);

qpast13 = DGregion_past(:,rx-1,ry+0,rz);
qpast14 = DGregion_past(:,rx+0,ry+0,rz);
qpast15 = DGregion_past(:,rx+1,ry+0,rz);

qpast16 = DGregion_past(:,rx-1,ry+1,rz);
qpast17 = DGregion_past(:,rx+0,ry+1,rz);
qpast18 = DGregion_past(:,rx+1,ry+1,rz);

qpast19 = DGregion_past(:,rx-1,ry-1,rz+1);
qpast20 = DGregion_past(:,rx+0,ry-1,rz+1);
qpast21 = DGregion_past(:,rx+1,ry-1,rz+1);

qpast22 = DGregion_past(:,rx-1,ry+0,rz+1);
qpast23 = DGregion_past(:,rx+0,ry+0,rz+1);
qpast24 = DGregion_past(:,rx+1,ry+0,rz+1);

qpast25 = DGregion_past(:,rx-1,ry+1,rz+1);
qpast26 = DGregion_past(:,rx+0,ry+1,rz+1);
qpast27 = DGregion_past(:,rx+1,ry+1,rz+1);



periodic = @(i,N) mod(i-1,N)+1; 
ivx = cellindex_center(1);
ivy = cellindex_center(2);
ivz = cellindex_center(3);

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

% PRESSURE
P01 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
P02 = auxiliary.total_pressure(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
P03 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));

P04 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz-1,Nv3));
P05 = auxiliary.total_pressure(:,ivx+0,ivy+0,periodic(ivz-1,Nv3));
P06 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz-1,Nv3));

P07 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
P08 = auxiliary.total_pressure(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
P09 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));

P10 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),ivz);
P11 = auxiliary.total_pressure(:,ivx+0,periodic(ivy-1,Nv2),ivz);
P12 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),ivz);

P13 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),ivy+0,ivz);
P14 = auxiliary.total_pressure(:,ivx+0,ivy+0,ivz);
P15 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),ivy+0,ivz);

P16 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),ivz);
P17 = auxiliary.total_pressure(:,ivx+0,periodic(ivy+1,Nv2),ivz);
P18 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),ivz);

P19 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
P20 = auxiliary.total_pressure(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
P21 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));

P22 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz+1,Nv3));
P23 = auxiliary.total_pressure(:,ivx+0,ivy+0,periodic(ivz+1,Nv3));
P24 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz+1,Nv3));

P25 = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
P26 = auxiliary.total_pressure(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
P27 = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));


DPDN01 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
DPDN02 = auxiliary.dpdn(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
DPDN03 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));

DPDN04 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz-1,Nv3));
DPDN05 = auxiliary.dpdn(:,ivx+0,ivy+0,periodic(ivz-1,Nv3));
DPDN06 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz-1,Nv3));

DPDN07 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
DPDN08 = auxiliary.dpdn(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
DPDN09 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));

DPDN10 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),ivz);
DPDN11 = auxiliary.dpdn(:,ivx+0,periodic(ivy-1,Nv2),ivz);
DPDN12 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),ivz);

DPDN13 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),ivy+0,ivz);
DPDN14 = auxiliary.dpdn(:,ivx+0,ivy+0,ivz);
DPDN15 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),ivy+0,ivz);

DPDN16 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),ivz);
DPDN17 = auxiliary.dpdn(:,ivx+0,periodic(ivy+1,Nv2),ivz);
DPDN18 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),ivz);

DPDN19 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
DPDN20 = auxiliary.dpdn(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
DPDN21 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));

DPDN22 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz+1,Nv3));
DPDN23 = auxiliary.dpdn(:,ivx+0,ivy+0,periodic(ivz+1,Nv3));
DPDN24 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz+1,Nv3));

DPDN25 = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
DPDN26 = auxiliary.dpdn(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
DPDN27 = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));


total_viscocity01 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
total_viscocity02 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));
total_viscocity03 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz-1,Nv3));

total_viscocity04 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz-1,Nv3));
total_viscocity05 = auxiliary.total_viscocity(:,ivx+0,ivy+0,periodic(ivz-1,Nv3));
total_viscocity06 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz-1,Nv3));

total_viscocity07 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
total_viscocity08 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));
total_viscocity09 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz-1,Nv3));

total_viscocity10 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),ivz);
total_viscocity11 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy-1,Nv2),ivz);
total_viscocity12 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),ivz);

total_viscocity13 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),ivy+0,ivz);
total_viscocity14 = auxiliary.total_viscocity(:,ivx+0,ivy+0,ivz);
total_viscocity15 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),ivy+0,ivz);

total_viscocity16 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),ivz);
total_viscocity17 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy+1,Nv2),ivz);
total_viscocity18 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),ivz);

total_viscocity19 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
total_viscocity20 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));
total_viscocity21 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy-1,Nv2),periodic(ivz+1,Nv3));

total_viscocity22 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),ivy+0,periodic(ivz+1,Nv3));
total_viscocity23 = auxiliary.total_viscocity(:,ivx+0,ivy+0,periodic(ivz+1,Nv3));
total_viscocity24 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),ivy+0,periodic(ivz+1,Nv3));

total_viscocity25 = auxiliary.total_viscocity(:,periodic(ivx-1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
total_viscocity26 = auxiliary.total_viscocity(:,ivx+0,periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));
total_viscocity27 = auxiliary.total_viscocity(:,periodic(ivx+1,Nv1),periodic(ivy+1,Nv2),periodic(ivz+1,Nv3));


[Jac_cell_111,residual_cell_111] = problem_jac_res_init(q01,qpast01,P01,DPDN01,total_viscocity01,auxiliary,data);
[Jac_cell_211,residual_cell_211] = problem_jac_res_init(q02,qpast02,P02,DPDN02,total_viscocity02,auxiliary,data);
[Jac_cell_311,residual_cell_311] = problem_jac_res_init(q03,qpast03,P03,DPDN03,total_viscocity03,auxiliary,data);
[Jac_cell_121,residual_cell_121] = problem_jac_res_init(q04,qpast04,P04,DPDN04,total_viscocity04,auxiliary,data);
[Jac_cell_221,residual_cell_221] = problem_jac_res_init(q05,qpast05,P05,DPDN05,total_viscocity05,auxiliary,data);
[Jac_cell_321,residual_cell_321] = problem_jac_res_init(q06,qpast06,P06,DPDN06,total_viscocity06,auxiliary,data);
[Jac_cell_131,residual_cell_131] = problem_jac_res_init(q07,qpast07,P07,DPDN07,total_viscocity07,auxiliary,data);
[Jac_cell_231,residual_cell_231] = problem_jac_res_init(q08,qpast08,P08,DPDN08,total_viscocity08,auxiliary,data);
[Jac_cell_331,residual_cell_331] = problem_jac_res_init(q09,qpast09,P09,DPDN09,total_viscocity09,auxiliary,data);

[Jac_cell_112,residual_cell_112] = problem_jac_res_init(q10,qpast10,P10,DPDN10,total_viscocity10,auxiliary,data);
[Jac_cell_212,residual_cell_212] = problem_jac_res_init(q11,qpast11,P11,DPDN11,total_viscocity11,auxiliary,data);
[Jac_cell_312,residual_cell_312] = problem_jac_res_init(q12,qpast12,P12,DPDN12,total_viscocity12,auxiliary,data);
[Jac_cell_122,residual_cell_122] = problem_jac_res_init(q13,qpast13,P13,DPDN13,total_viscocity13,auxiliary,data);
[Jac_cell_222,residual_cell_222] = problem_jac_res_init(q14,qpast14,P14,DPDN14,total_viscocity14,auxiliary,data);
[Jac_cell_322,residual_cell_322] = problem_jac_res_init(q15,qpast15,P15,DPDN15,total_viscocity15,auxiliary,data);
[Jac_cell_132,residual_cell_132] = problem_jac_res_init(q16,qpast16,P16,DPDN16,total_viscocity16,auxiliary,data);
[Jac_cell_232,residual_cell_232] = problem_jac_res_init(q17,qpast17,P17,DPDN17,total_viscocity17,auxiliary,data);
[Jac_cell_332,residual_cell_332] = problem_jac_res_init(q18,qpast18,P18,DPDN18,total_viscocity18,auxiliary,data);

[Jac_cell_113,residual_cell_113] = problem_jac_res_init(q19,qpast19,P19,DPDN19,total_viscocity19,auxiliary,data);
[Jac_cell_213,residual_cell_213] = problem_jac_res_init(q20,qpast20,P20,DPDN20,total_viscocity20,auxiliary,data);
[Jac_cell_313,residual_cell_313] = problem_jac_res_init(q21,qpast21,P21,DPDN21,total_viscocity21,auxiliary,data);
[Jac_cell_123,residual_cell_123] = problem_jac_res_init(q22,qpast22,P22,DPDN22,total_viscocity22,auxiliary,data);
[Jac_cell_223,residual_cell_223] = problem_jac_res_init(q23,qpast23,P23,DPDN23,total_viscocity23,auxiliary,data);
[Jac_cell_323,residual_cell_323] = problem_jac_res_init(q24,qpast24,P24,DPDN24,total_viscocity24,auxiliary,data);
[Jac_cell_133,residual_cell_133] = problem_jac_res_init(q25,qpast25,P25,DPDN25,total_viscocity25,auxiliary,data);
[Jac_cell_233,residual_cell_233] = problem_jac_res_init(q26,qpast26,P26,DPDN26,total_viscocity26,auxiliary,data);
[Jac_cell_333,residual_cell_333] = problem_jac_res_init(q27,qpast27,P27,DPDN27,total_viscocity27,auxiliary,data);


%RESIDUAL
east_flux_111= 0 ; %01
east_flux_211= 0 ; %02
east_flux_311= 0 ; %03
east_flux_121= 0 ; %04
east_flux_221= 0 ; %05
east_flux_321= 0 ; %06
east_flux_131= 0 ; %07
east_flux_231= 0 ; %08
east_flux_331= 0 ; %09
east_flux_112= 0 ; %10
east_flux_212= 0 ; %11
east_flux_312= 0 ; %12
east_flux_122= 0 ; %13
east_flux_222= 0 ; %14
east_flux_322= 0 ; %15
east_flux_132= 0 ; %16
east_flux_232= 0 ; %17
east_flux_332= 0 ; %18
east_flux_113= 0 ; %19
east_flux_213= 0 ; %20
east_flux_313= 0 ; %21
east_flux_123= 0 ; %22
east_flux_223= 0 ; %23
east_flux_323= 0 ; %24
east_flux_133= 0 ; %25
east_flux_233= 0 ; %26
east_flux_333= 0 ; %27
west_flux_111= 0 ; %01
west_flux_211= 0 ; %02
west_flux_311= 0 ; %03
west_flux_121= 0 ; %04
west_flux_221= 0 ; %05
west_flux_321= 0 ; %06
west_flux_131= 0 ; %07
west_flux_231= 0 ; %08
west_flux_331= 0 ; %09
west_flux_112= 0 ; %10
west_flux_212= 0 ; %11
west_flux_312= 0 ; %12
west_flux_122= 0 ; %13
west_flux_222= 0 ; %14
west_flux_322= 0 ; %15
west_flux_132= 0 ; %16
west_flux_232= 0 ; %17
west_flux_332= 0 ; %18
west_flux_113= 0 ; %19
west_flux_213= 0 ; %20
west_flux_313= 0 ; %21
west_flux_123= 0 ; %22
west_flux_223= 0 ; %23
west_flux_323= 0 ; %24
west_flux_133= 0 ; %25
west_flux_233= 0 ; %26
west_flux_333= 0 ; %27
nort_flux_111= 0 ; %01
nort_flux_211= 0 ; %02
nort_flux_311= 0 ; %03
nort_flux_121= 0 ; %04
nort_flux_221= 0 ; %05
nort_flux_321= 0 ; %06
nort_flux_131= 0 ; %07
nort_flux_231= 0 ; %08
nort_flux_331= 0 ; %09
nort_flux_112= 0 ; %10
nort_flux_212= 0 ; %11
nort_flux_312= 0 ; %12
nort_flux_122= 0 ; %13
nort_flux_222= 0 ; %14
nort_flux_322= 0 ; %15
nort_flux_132= 0 ; %16
nort_flux_232= 0 ; %17
nort_flux_332= 0 ; %18
nort_flux_113= 0 ; %19
nort_flux_213= 0 ; %20
nort_flux_313= 0 ; %21
nort_flux_123= 0 ; %22
nort_flux_223= 0 ; %23
nort_flux_323= 0 ; %24
nort_flux_133= 0 ; %25
nort_flux_233= 0 ; %26
nort_flux_333= 0 ; %27
sout_flux_111= 0 ; %01
sout_flux_211= 0 ; %02
sout_flux_311= 0 ; %03
sout_flux_121= 0 ; %04
sout_flux_221= 0 ; %05
sout_flux_321= 0 ; %06
sout_flux_131= 0 ; %07
sout_flux_231= 0 ; %08
sout_flux_331= 0 ; %09
sout_flux_112= 0 ; %10
sout_flux_212= 0 ; %11
sout_flux_312= 0 ; %12
sout_flux_122= 0 ; %13
sout_flux_222= 0 ; %14
sout_flux_322= 0 ; %15
sout_flux_132= 0 ; %16
sout_flux_232= 0 ; %17
sout_flux_332= 0 ; %18
sout_flux_113= 0 ; %19
sout_flux_213= 0 ; %20
sout_flux_313= 0 ; %21
sout_flux_123= 0 ; %22
sout_flux_223= 0 ; %23
sout_flux_323= 0 ; %24
sout_flux_133= 0 ; %25
sout_flux_233= 0 ; %26
sout_flux_333= 0 ; %27
uppr_flux_111= 0 ; %01
uppr_flux_211= 0 ; %02
uppr_flux_311= 0 ; %03
uppr_flux_121= 0 ; %04
uppr_flux_221= 0 ; %05
uppr_flux_321= 0 ; %06
uppr_flux_131= 0 ; %07
uppr_flux_231= 0 ; %08
uppr_flux_331= 0 ; %09
uppr_flux_112= 0 ; %10
uppr_flux_212= 0 ; %11
uppr_flux_312= 0 ; %12
uppr_flux_122= 0 ; %13
uppr_flux_222= 0 ; %14
uppr_flux_322= 0 ; %15
uppr_flux_132= 0 ; %16
uppr_flux_232= 0 ; %17
uppr_flux_332= 0 ; %18
uppr_flux_113= 0 ; %19
uppr_flux_213= 0 ; %20
uppr_flux_313= 0 ; %21
uppr_flux_123= 0 ; %22
uppr_flux_223= 0 ; %23
uppr_flux_323= 0 ; %24
uppr_flux_133= 0 ; %25
uppr_flux_233= 0 ; %26
uppr_flux_333= 0 ; %27
down_flux_111= 0 ; %01
down_flux_211= 0 ; %02
down_flux_311= 0 ; %03
down_flux_121= 0 ; %04
down_flux_221= 0 ; %05
down_flux_321= 0 ; %06
down_flux_131= 0 ; %07
down_flux_231= 0 ; %08
down_flux_331= 0 ; %09
down_flux_112= 0 ; %10
down_flux_212= 0 ; %11
down_flux_312= 0 ; %12
down_flux_122= 0 ; %13
down_flux_222= 0 ; %14
down_flux_322= 0 ; %15
down_flux_132= 0 ; %16
down_flux_232= 0 ; %17
down_flux_332= 0 ; %18
down_flux_113= 0 ; %19
down_flux_213= 0 ; %20
down_flux_313= 0 ; %21
down_flux_123= 0 ; %22
down_flux_223= 0 ; %23
down_flux_323= 0 ; %24
down_flux_133= 0 ; %25
down_flux_233= 0 ; %26
down_flux_333= 0 ; %27

east_trunc_111= 0 ; %01
east_trunc_211= 0 ; %02
east_trunc_311= 0 ; %03
east_trunc_121= 0 ; %04
east_trunc_221= 0 ; %05
east_trunc_321= 0 ; %06
east_trunc_131= 0 ; %07
east_trunc_231= 0 ; %08
east_trunc_331= 0 ; %09
east_trunc_112= 0 ; %10
east_trunc_212= 0 ; %11
east_trunc_312= 0 ; %12
east_trunc_122= 0 ; %13
east_trunc_222= 0 ; %14
east_trunc_322= 0 ; %15
east_trunc_132= 0 ; %16
east_trunc_232= 0 ; %17
east_trunc_332= 0 ; %18
east_trunc_113= 0 ; %19
east_trunc_213= 0 ; %20
east_trunc_313= 0 ; %21
east_trunc_123= 0 ; %22
east_trunc_223= 0 ; %23
east_trunc_323= 0 ; %24
east_trunc_133= 0 ; %25
east_trunc_233= 0 ; %26
east_trunc_333= 0 ; %27
west_trunc_111= 0 ; %01
west_trunc_211= 0 ; %02
west_trunc_311= 0 ; %03
west_trunc_121= 0 ; %04
west_trunc_221= 0 ; %05
west_trunc_321= 0 ; %06
west_trunc_131= 0 ; %07
west_trunc_231= 0 ; %08
west_trunc_331= 0 ; %09
west_trunc_112= 0 ; %10
west_trunc_212= 0 ; %11
west_trunc_312= 0 ; %12
west_trunc_122= 0 ; %13
west_trunc_222= 0 ; %14
west_trunc_322= 0 ; %15
west_trunc_132= 0 ; %16
west_trunc_232= 0 ; %17
west_trunc_332= 0 ; %18
west_trunc_113= 0 ; %19
west_trunc_213= 0 ; %20
west_trunc_313= 0 ; %21
west_trunc_123= 0 ; %22
west_trunc_223= 0 ; %23
west_trunc_323= 0 ; %24
west_trunc_133= 0 ; %25
west_trunc_233= 0 ; %26
west_trunc_333= 0 ; %27
nort_trunc_111= 0 ; %01
nort_trunc_211= 0 ; %02
nort_trunc_311= 0 ; %03
nort_trunc_121= 0 ; %04
nort_trunc_221= 0 ; %05
nort_trunc_321= 0 ; %06
nort_trunc_131= 0 ; %07
nort_trunc_231= 0 ; %08
nort_trunc_331= 0 ; %09
nort_trunc_112= 0 ; %10
nort_trunc_212= 0 ; %11
nort_trunc_312= 0 ; %12
nort_trunc_122= 0 ; %13
nort_trunc_222= 0 ; %14
nort_trunc_322= 0 ; %15
nort_trunc_132= 0 ; %16
nort_trunc_232= 0 ; %17
nort_trunc_332= 0 ; %18
nort_trunc_113= 0 ; %19
nort_trunc_213= 0 ; %20
nort_trunc_313= 0 ; %21
nort_trunc_123= 0 ; %22
nort_trunc_223= 0 ; %23
nort_trunc_323= 0 ; %24
nort_trunc_133= 0 ; %25
nort_trunc_233= 0 ; %26
nort_trunc_333= 0 ; %27
sout_trunc_111= 0 ; %01
sout_trunc_211= 0 ; %02
sout_trunc_311= 0 ; %03
sout_trunc_121= 0 ; %04
sout_trunc_221= 0 ; %05
sout_trunc_321= 0 ; %06
sout_trunc_131= 0 ; %07
sout_trunc_231= 0 ; %08
sout_trunc_331= 0 ; %09
sout_trunc_112= 0 ; %10
sout_trunc_212= 0 ; %11
sout_trunc_312= 0 ; %12
sout_trunc_122= 0 ; %13
sout_trunc_222= 0 ; %14
sout_trunc_322= 0 ; %15
sout_trunc_132= 0 ; %16
sout_trunc_232= 0 ; %17
sout_trunc_332= 0 ; %18
sout_trunc_113= 0 ; %19
sout_trunc_213= 0 ; %20
sout_trunc_313= 0 ; %21
sout_trunc_123= 0 ; %22
sout_trunc_223= 0 ; %23
sout_trunc_323= 0 ; %24
sout_trunc_133= 0 ; %25
sout_trunc_233= 0 ; %26
sout_trunc_333= 0 ; %27
uppr_trunc_111= 0 ; %01
uppr_trunc_211= 0 ; %02
uppr_trunc_311= 0 ; %03
uppr_trunc_121= 0 ; %04
uppr_trunc_221= 0 ; %05
uppr_trunc_321= 0 ; %06
uppr_trunc_131= 0 ; %07
uppr_trunc_231= 0 ; %08
uppr_trunc_331= 0 ; %09
uppr_trunc_112= 0 ; %10
uppr_trunc_212= 0 ; %11
uppr_trunc_312= 0 ; %12
uppr_trunc_122= 0 ; %13
uppr_trunc_222= 0 ; %14
uppr_trunc_322= 0 ; %15
uppr_trunc_132= 0 ; %16
uppr_trunc_232= 0 ; %17
uppr_trunc_332= 0 ; %18
uppr_trunc_113= 0 ; %19
uppr_trunc_213= 0 ; %20
uppr_trunc_313= 0 ; %21
uppr_trunc_123= 0 ; %22
uppr_trunc_223= 0 ; %23
uppr_trunc_323= 0 ; %24
uppr_trunc_133= 0 ; %25
uppr_trunc_233= 0 ; %26
uppr_trunc_333= 0 ; %27
down_trunc_111= 0 ; %01
down_trunc_211= 0 ; %02
down_trunc_311= 0 ; %03
down_trunc_121= 0 ; %04
down_trunc_221= 0 ; %05
down_trunc_321= 0 ; %06
down_trunc_131= 0 ; %07
down_trunc_231= 0 ; %08
down_trunc_331= 0 ; %09
down_trunc_112= 0 ; %10
down_trunc_212= 0 ; %11
down_trunc_312= 0 ; %12
down_trunc_122= 0 ; %13
down_trunc_222= 0 ; %14
down_trunc_322= 0 ; %15
down_trunc_132= 0 ; %16
down_trunc_232= 0 ; %17
down_trunc_332= 0 ; %18
down_trunc_113= 0 ; %19
down_trunc_213= 0 ; %20
down_trunc_313= 0 ; %21
down_trunc_123= 0 ; %22
down_trunc_223= 0 ; %23
down_trunc_323= 0 ; %24
down_trunc_133= 0 ; %25
down_trunc_233= 0 ; %26
down_trunc_333= 0 ; %27


% JACOBIAN 
Jac_east_flux_111= 0 ; %01
Jac_east_flux_211= 0 ; %02
Jac_east_flux_311= 0 ; %03
Jac_east_flux_121= 0 ; %04
Jac_east_flux_221= 0 ; %05
Jac_east_flux_321= 0 ; %06
Jac_east_flux_131= 0 ; %07
Jac_east_flux_231= 0 ; %08
Jac_east_flux_331= 0 ; %09
Jac_east_flux_112= 0 ; %10
Jac_east_flux_212= 0 ; %11
Jac_east_flux_312= 0 ; %12
Jac_east_flux_122= 0 ; %13
Jac_east_flux_222= 0 ; %14
Jac_east_flux_322= 0 ; %15
Jac_east_flux_132= 0 ; %16
Jac_east_flux_232= 0 ; %17
Jac_east_flux_332= 0 ; %18
Jac_east_flux_113= 0 ; %19
Jac_east_flux_213= 0 ; %20
Jac_east_flux_313= 0 ; %21
Jac_east_flux_123= 0 ; %22
Jac_east_flux_223= 0 ; %23
Jac_east_flux_323= 0 ; %24
Jac_east_flux_133= 0 ; %25
Jac_east_flux_233= 0 ; %26
Jac_east_flux_333= 0 ; %27
Jac_west_flux_111= 0 ; %01
Jac_west_flux_211= 0 ; %02
Jac_west_flux_311= 0 ; %03
Jac_west_flux_121= 0 ; %04
Jac_west_flux_221= 0 ; %05
Jac_west_flux_321= 0 ; %06
Jac_west_flux_131= 0 ; %07
Jac_west_flux_231= 0 ; %08
Jac_west_flux_331= 0 ; %09
Jac_west_flux_112= 0 ; %10
Jac_west_flux_212= 0 ; %11
Jac_west_flux_312= 0 ; %12
Jac_west_flux_122= 0 ; %13
Jac_west_flux_222= 0 ; %14
Jac_west_flux_322= 0 ; %15
Jac_west_flux_132= 0 ; %16
Jac_west_flux_232= 0 ; %17
Jac_west_flux_332= 0 ; %18
Jac_west_flux_113= 0 ; %19
Jac_west_flux_213= 0 ; %20
Jac_west_flux_313= 0 ; %21
Jac_west_flux_123= 0 ; %22
Jac_west_flux_223= 0 ; %23
Jac_west_flux_323= 0 ; %24
Jac_west_flux_133= 0 ; %25
Jac_west_flux_233= 0 ; %26
Jac_west_flux_333= 0 ; %27
Jac_nort_flux_111= 0 ; %01
Jac_nort_flux_211= 0 ; %02
Jac_nort_flux_311= 0 ; %03
Jac_nort_flux_121= 0 ; %04
Jac_nort_flux_221= 0 ; %05
Jac_nort_flux_321= 0 ; %06
Jac_nort_flux_131= 0 ; %07
Jac_nort_flux_231= 0 ; %08
Jac_nort_flux_331= 0 ; %09
Jac_nort_flux_112= 0 ; %10
Jac_nort_flux_212= 0 ; %11
Jac_nort_flux_312= 0 ; %12
Jac_nort_flux_122= 0 ; %13
Jac_nort_flux_222= 0 ; %14
Jac_nort_flux_322= 0 ; %15
Jac_nort_flux_132= 0 ; %16
Jac_nort_flux_232= 0 ; %17
Jac_nort_flux_332= 0 ; %18
Jac_nort_flux_113= 0 ; %19
Jac_nort_flux_213= 0 ; %20
Jac_nort_flux_313= 0 ; %21
Jac_nort_flux_123= 0 ; %22
Jac_nort_flux_223= 0 ; %23
Jac_nort_flux_323= 0 ; %24
Jac_nort_flux_133= 0 ; %25
Jac_nort_flux_233= 0 ; %26
Jac_nort_flux_333= 0 ; %27
Jac_sout_flux_111= 0 ; %01
Jac_sout_flux_211= 0 ; %02
Jac_sout_flux_311= 0 ; %03
Jac_sout_flux_121= 0 ; %04
Jac_sout_flux_221= 0 ; %05
Jac_sout_flux_321= 0 ; %06
Jac_sout_flux_131= 0 ; %07
Jac_sout_flux_231= 0 ; %08
Jac_sout_flux_331= 0 ; %09
Jac_sout_flux_112= 0 ; %10
Jac_sout_flux_212= 0 ; %11
Jac_sout_flux_312= 0 ; %12
Jac_sout_flux_122= 0 ; %13
Jac_sout_flux_222= 0 ; %14
Jac_sout_flux_322= 0 ; %15
Jac_sout_flux_132= 0 ; %16
Jac_sout_flux_232= 0 ; %17
Jac_sout_flux_332= 0 ; %18
Jac_sout_flux_113= 0 ; %19
Jac_sout_flux_213= 0 ; %20
Jac_sout_flux_313= 0 ; %21
Jac_sout_flux_123= 0 ; %22
Jac_sout_flux_223= 0 ; %23
Jac_sout_flux_323= 0 ; %24
Jac_sout_flux_133= 0 ; %25
Jac_sout_flux_233= 0 ; %26
Jac_sout_flux_333= 0 ; %27
Jac_uppr_flux_111= 0 ; %01
Jac_uppr_flux_211= 0 ; %02
Jac_uppr_flux_311= 0 ; %03
Jac_uppr_flux_121= 0 ; %04
Jac_uppr_flux_221= 0 ; %05
Jac_uppr_flux_321= 0 ; %06
Jac_uppr_flux_131= 0 ; %07
Jac_uppr_flux_231= 0 ; %08
Jac_uppr_flux_331= 0 ; %09
Jac_uppr_flux_112= 0 ; %10
Jac_uppr_flux_212= 0 ; %11
Jac_uppr_flux_312= 0 ; %12
Jac_uppr_flux_122= 0 ; %13
Jac_uppr_flux_222= 0 ; %14
Jac_uppr_flux_322= 0 ; %15
Jac_uppr_flux_132= 0 ; %16
Jac_uppr_flux_232= 0 ; %17
Jac_uppr_flux_332= 0 ; %18
Jac_uppr_flux_113= 0 ; %19
Jac_uppr_flux_213= 0 ; %20
Jac_uppr_flux_313= 0 ; %21
Jac_uppr_flux_123= 0 ; %22
Jac_uppr_flux_223= 0 ; %23
Jac_uppr_flux_323= 0 ; %24
Jac_uppr_flux_133= 0 ; %25
Jac_uppr_flux_233= 0 ; %26
Jac_uppr_flux_333= 0 ; %27
Jac_down_flux_111= 0 ; %01
Jac_down_flux_211= 0 ; %02
Jac_down_flux_311= 0 ; %03
Jac_down_flux_121= 0 ; %04
Jac_down_flux_221= 0 ; %05
Jac_down_flux_321= 0 ; %06
Jac_down_flux_131= 0 ; %07
Jac_down_flux_231= 0 ; %08
Jac_down_flux_331= 0 ; %09
Jac_down_flux_112= 0 ; %10
Jac_down_flux_212= 0 ; %11
Jac_down_flux_312= 0 ; %12
Jac_down_flux_122= 0 ; %13
Jac_down_flux_222= 0 ; %14
Jac_down_flux_322= 0 ; %15
Jac_down_flux_132= 0 ; %16
Jac_down_flux_232= 0 ; %17
Jac_down_flux_332= 0 ; %18
Jac_down_flux_113= 0 ; %19
Jac_down_flux_213= 0 ; %20
Jac_down_flux_313= 0 ; %21
Jac_down_flux_123= 0 ; %22
Jac_down_flux_223= 0 ; %23
Jac_down_flux_323= 0 ; %24
Jac_down_flux_133= 0 ; %25
Jac_down_flux_233= 0 ; %26
Jac_down_flux_333= 0 ; %27

Jac_east_othr_111= 0 ; %01
Jac_east_othr_211= 0 ; %02
Jac_east_othr_311= 0 ; %03
Jac_east_othr_121= 0 ; %04
Jac_east_othr_221= 0 ; %05
Jac_east_othr_321= 0 ; %06
Jac_east_othr_131= 0 ; %07
Jac_east_othr_231= 0 ; %08
Jac_east_othr_331= 0 ; %09
Jac_east_othr_112= 0 ; %10
Jac_east_othr_212= 0 ; %11
Jac_east_othr_312= 0 ; %12
Jac_east_othr_122= 0 ; %13
Jac_east_othr_222= 0 ; %14
Jac_east_othr_322= 0 ; %15
Jac_east_othr_132= 0 ; %16
Jac_east_othr_232= 0 ; %17
Jac_east_othr_332= 0 ; %18
Jac_east_othr_113= 0 ; %19
Jac_east_othr_213= 0 ; %20
Jac_east_othr_313= 0 ; %21
Jac_east_othr_123= 0 ; %22
Jac_east_othr_223= 0 ; %23
Jac_east_othr_323= 0 ; %24
Jac_east_othr_133= 0 ; %25
Jac_east_othr_233= 0 ; %26
Jac_east_othr_333= 0 ; %27
Jac_west_othr_111= 0 ; %01
Jac_west_othr_211= 0 ; %02
Jac_west_othr_311= 0 ; %03
Jac_west_othr_121= 0 ; %04
Jac_west_othr_221= 0 ; %05
Jac_west_othr_321= 0 ; %06
Jac_west_othr_131= 0 ; %07
Jac_west_othr_231= 0 ; %08
Jac_west_othr_331= 0 ; %09
Jac_west_othr_112= 0 ; %10
Jac_west_othr_212= 0 ; %11
Jac_west_othr_312= 0 ; %12
Jac_west_othr_122= 0 ; %13
Jac_west_othr_222= 0 ; %14
Jac_west_othr_322= 0 ; %15
Jac_west_othr_132= 0 ; %16
Jac_west_othr_232= 0 ; %17
Jac_west_othr_332= 0 ; %18
Jac_west_othr_113= 0 ; %19
Jac_west_othr_213= 0 ; %20
Jac_west_othr_313= 0 ; %21
Jac_west_othr_123= 0 ; %22
Jac_west_othr_223= 0 ; %23
Jac_west_othr_323= 0 ; %24
Jac_west_othr_133= 0 ; %25
Jac_west_othr_233= 0 ; %26
Jac_west_othr_333= 0 ; %27
Jac_nort_othr_111= 0 ; %01
Jac_nort_othr_211= 0 ; %02
Jac_nort_othr_311= 0 ; %03
Jac_nort_othr_121= 0 ; %04
Jac_nort_othr_221= 0 ; %05
Jac_nort_othr_321= 0 ; %06
Jac_nort_othr_131= 0 ; %07
Jac_nort_othr_231= 0 ; %08
Jac_nort_othr_331= 0 ; %09
Jac_nort_othr_112= 0 ; %10
Jac_nort_othr_212= 0 ; %11
Jac_nort_othr_312= 0 ; %12
Jac_nort_othr_122= 0 ; %13
Jac_nort_othr_222= 0 ; %14
Jac_nort_othr_322= 0 ; %15
Jac_nort_othr_132= 0 ; %16
Jac_nort_othr_232= 0 ; %17
Jac_nort_othr_332= 0 ; %18
Jac_nort_othr_113= 0 ; %19
Jac_nort_othr_213= 0 ; %20
Jac_nort_othr_313= 0 ; %21
Jac_nort_othr_123= 0 ; %22
Jac_nort_othr_223= 0 ; %23
Jac_nort_othr_323= 0 ; %24
Jac_nort_othr_133= 0 ; %25
Jac_nort_othr_233= 0 ; %26
Jac_nort_othr_333= 0 ; %27
Jac_sout_othr_111= 0 ; %01
Jac_sout_othr_211= 0 ; %02
Jac_sout_othr_311= 0 ; %03
Jac_sout_othr_121= 0 ; %04
Jac_sout_othr_221= 0 ; %05
Jac_sout_othr_321= 0 ; %06
Jac_sout_othr_131= 0 ; %07
Jac_sout_othr_231= 0 ; %08
Jac_sout_othr_331= 0 ; %09
Jac_sout_othr_112= 0 ; %10
Jac_sout_othr_212= 0 ; %11
Jac_sout_othr_312= 0 ; %12
Jac_sout_othr_122= 0 ; %13
Jac_sout_othr_222= 0 ; %14
Jac_sout_othr_322= 0 ; %15
Jac_sout_othr_132= 0 ; %16
Jac_sout_othr_232= 0 ; %17
Jac_sout_othr_332= 0 ; %18
Jac_sout_othr_113= 0 ; %19
Jac_sout_othr_213= 0 ; %20
Jac_sout_othr_313= 0 ; %21
Jac_sout_othr_123= 0 ; %22
Jac_sout_othr_223= 0 ; %23
Jac_sout_othr_323= 0 ; %24
Jac_sout_othr_133= 0 ; %25
Jac_sout_othr_233= 0 ; %26
Jac_sout_othr_333= 0 ; %27
Jac_uppr_othr_111= 0 ; %01
Jac_uppr_othr_211= 0 ; %02
Jac_uppr_othr_311= 0 ; %03
Jac_uppr_othr_121= 0 ; %04
Jac_uppr_othr_221= 0 ; %05
Jac_uppr_othr_321= 0 ; %06
Jac_uppr_othr_131= 0 ; %07
Jac_uppr_othr_231= 0 ; %08
Jac_uppr_othr_331= 0 ; %09
Jac_uppr_othr_112= 0 ; %10
Jac_uppr_othr_212= 0 ; %11
Jac_uppr_othr_312= 0 ; %12
Jac_uppr_othr_122= 0 ; %13
Jac_uppr_othr_222= 0 ; %14
Jac_uppr_othr_322= 0 ; %15
Jac_uppr_othr_132= 0 ; %16
Jac_uppr_othr_232= 0 ; %17
Jac_uppr_othr_332= 0 ; %18
Jac_uppr_othr_113= 0 ; %19
Jac_uppr_othr_213= 0 ; %20
Jac_uppr_othr_313= 0 ; %21
Jac_uppr_othr_123= 0 ; %22
Jac_uppr_othr_223= 0 ; %23
Jac_uppr_othr_323= 0 ; %24
Jac_uppr_othr_133= 0 ; %25
Jac_uppr_othr_233= 0 ; %26
Jac_uppr_othr_333= 0 ; %27
Jac_down_othr_111= 0 ; %01
Jac_down_othr_211= 0 ; %02
Jac_down_othr_311= 0 ; %03
Jac_down_othr_121= 0 ; %04
Jac_down_othr_221= 0 ; %05
Jac_down_othr_321= 0 ; %06
Jac_down_othr_131= 0 ; %07
Jac_down_othr_231= 0 ; %08
Jac_down_othr_331= 0 ; %09
Jac_down_othr_112= 0 ; %10
Jac_down_othr_212= 0 ; %11
Jac_down_othr_312= 0 ; %12
Jac_down_othr_122= 0 ; %13
Jac_down_othr_222= 0 ; %14
Jac_down_othr_322= 0 ; %15
Jac_down_othr_132= 0 ; %16
Jac_down_othr_232= 0 ; %17
Jac_down_othr_332= 0 ; %18
Jac_down_othr_113= 0 ; %19
Jac_down_othr_213= 0 ; %20
Jac_down_othr_313= 0 ; %21
Jac_down_othr_123= 0 ; %22
Jac_down_othr_223= 0 ; %23
Jac_down_othr_323= 0 ; %24
Jac_down_othr_133= 0 ; %25
Jac_down_othr_233= 0 ; %26
Jac_down_othr_333= 0 ; %27

Jac_east_trunc_111= 0 ; %01
Jac_east_trunc_211= 0 ; %02
Jac_east_trunc_311= 0 ; %03
Jac_east_trunc_121= 0 ; %04
Jac_east_trunc_221= 0 ; %05
Jac_east_trunc_321= 0 ; %06
Jac_east_trunc_131= 0 ; %07
Jac_east_trunc_231= 0 ; %08
Jac_east_trunc_331= 0 ; %09
Jac_east_trunc_112= 0 ; %10
Jac_east_trunc_212= 0 ; %11
Jac_east_trunc_312= 0 ; %12
Jac_east_trunc_122= 0 ; %13
Jac_east_trunc_222= 0 ; %14
Jac_east_trunc_322= 0 ; %15
Jac_east_trunc_132= 0 ; %16
Jac_east_trunc_232= 0 ; %17
Jac_east_trunc_332= 0 ; %18
Jac_east_trunc_113= 0 ; %19
Jac_east_trunc_213= 0 ; %20
Jac_east_trunc_313= 0 ; %21
Jac_east_trunc_123= 0 ; %22
Jac_east_trunc_223= 0 ; %23
Jac_east_trunc_323= 0 ; %24
Jac_east_trunc_133= 0 ; %25
Jac_east_trunc_233= 0 ; %26
Jac_east_trunc_333= 0 ; %27
Jac_west_trunc_111= 0 ; %01
Jac_west_trunc_211= 0 ; %02
Jac_west_trunc_311= 0 ; %03
Jac_west_trunc_121= 0 ; %04
Jac_west_trunc_221= 0 ; %05
Jac_west_trunc_321= 0 ; %06
Jac_west_trunc_131= 0 ; %07
Jac_west_trunc_231= 0 ; %08
Jac_west_trunc_331= 0 ; %09
Jac_west_trunc_112= 0 ; %10
Jac_west_trunc_212= 0 ; %11
Jac_west_trunc_312= 0 ; %12
Jac_west_trunc_122= 0 ; %13
Jac_west_trunc_222= 0 ; %14
Jac_west_trunc_322= 0 ; %15
Jac_west_trunc_132= 0 ; %16
Jac_west_trunc_232= 0 ; %17
Jac_west_trunc_332= 0 ; %18
Jac_west_trunc_113= 0 ; %19
Jac_west_trunc_213= 0 ; %20
Jac_west_trunc_313= 0 ; %21
Jac_west_trunc_123= 0 ; %22
Jac_west_trunc_223= 0 ; %23
Jac_west_trunc_323= 0 ; %24
Jac_west_trunc_133= 0 ; %25
Jac_west_trunc_233= 0 ; %26
Jac_west_trunc_333= 0 ; %27
Jac_nort_trunc_111= 0 ; %01
Jac_nort_trunc_211= 0 ; %02
Jac_nort_trunc_311= 0 ; %03
Jac_nort_trunc_121= 0 ; %04
Jac_nort_trunc_221= 0 ; %05
Jac_nort_trunc_321= 0 ; %06
Jac_nort_trunc_131= 0 ; %07
Jac_nort_trunc_231= 0 ; %08
Jac_nort_trunc_331= 0 ; %09
Jac_nort_trunc_112= 0 ; %10
Jac_nort_trunc_212= 0 ; %11
Jac_nort_trunc_312= 0 ; %12
Jac_nort_trunc_122= 0 ; %13
Jac_nort_trunc_222= 0 ; %14
Jac_nort_trunc_322= 0 ; %15
Jac_nort_trunc_132= 0 ; %16
Jac_nort_trunc_232= 0 ; %17
Jac_nort_trunc_332= 0 ; %18
Jac_nort_trunc_113= 0 ; %19
Jac_nort_trunc_213= 0 ; %20
Jac_nort_trunc_313= 0 ; %21
Jac_nort_trunc_123= 0 ; %22
Jac_nort_trunc_223= 0 ; %23
Jac_nort_trunc_323= 0 ; %24
Jac_nort_trunc_133= 0 ; %25
Jac_nort_trunc_233= 0 ; %26
Jac_nort_trunc_333= 0 ; %27
Jac_sout_trunc_111= 0 ; %01
Jac_sout_trunc_211= 0 ; %02
Jac_sout_trunc_311= 0 ; %03
Jac_sout_trunc_121= 0 ; %04
Jac_sout_trunc_221= 0 ; %05
Jac_sout_trunc_321= 0 ; %06
Jac_sout_trunc_131= 0 ; %07
Jac_sout_trunc_231= 0 ; %08
Jac_sout_trunc_331= 0 ; %09
Jac_sout_trunc_112= 0 ; %10
Jac_sout_trunc_212= 0 ; %11
Jac_sout_trunc_312= 0 ; %12
Jac_sout_trunc_122= 0 ; %13
Jac_sout_trunc_222= 0 ; %14
Jac_sout_trunc_322= 0 ; %15
Jac_sout_trunc_132= 0 ; %16
Jac_sout_trunc_232= 0 ; %17
Jac_sout_trunc_332= 0 ; %18
Jac_sout_trunc_113= 0 ; %19
Jac_sout_trunc_213= 0 ; %20
Jac_sout_trunc_313= 0 ; %21
Jac_sout_trunc_123= 0 ; %22
Jac_sout_trunc_223= 0 ; %23
Jac_sout_trunc_323= 0 ; %24
Jac_sout_trunc_133= 0 ; %25
Jac_sout_trunc_233= 0 ; %26
Jac_sout_trunc_333= 0 ; %27
Jac_uppr_trunc_111= 0 ; %01
Jac_uppr_trunc_211= 0 ; %02
Jac_uppr_trunc_311= 0 ; %03
Jac_uppr_trunc_121= 0 ; %04
Jac_uppr_trunc_221= 0 ; %05
Jac_uppr_trunc_321= 0 ; %06
Jac_uppr_trunc_131= 0 ; %07
Jac_uppr_trunc_231= 0 ; %08
Jac_uppr_trunc_331= 0 ; %09
Jac_uppr_trunc_112= 0 ; %10
Jac_uppr_trunc_212= 0 ; %11
Jac_uppr_trunc_312= 0 ; %12
Jac_uppr_trunc_122= 0 ; %13
Jac_uppr_trunc_222= 0 ; %14
Jac_uppr_trunc_322= 0 ; %15
Jac_uppr_trunc_132= 0 ; %16
Jac_uppr_trunc_232= 0 ; %17
Jac_uppr_trunc_332= 0 ; %18
Jac_uppr_trunc_113= 0 ; %19
Jac_uppr_trunc_213= 0 ; %20
Jac_uppr_trunc_313= 0 ; %21
Jac_uppr_trunc_123= 0 ; %22
Jac_uppr_trunc_223= 0 ; %23
Jac_uppr_trunc_323= 0 ; %24
Jac_uppr_trunc_133= 0 ; %25
Jac_uppr_trunc_233= 0 ; %26
Jac_uppr_trunc_333= 0 ; %27
Jac_down_trunc_111= 0 ; %01
Jac_down_trunc_211= 0 ; %02
Jac_down_trunc_311= 0 ; %03
Jac_down_trunc_121= 0 ; %04
Jac_down_trunc_221= 0 ; %05
Jac_down_trunc_321= 0 ; %06
Jac_down_trunc_131= 0 ; %07
Jac_down_trunc_231= 0 ; %08
Jac_down_trunc_331= 0 ; %09
Jac_down_trunc_112= 0 ; %10
Jac_down_trunc_212= 0 ; %11
Jac_down_trunc_312= 0 ; %12
Jac_down_trunc_122= 0 ; %13
Jac_down_trunc_222= 0 ; %14
Jac_down_trunc_322= 0 ; %15
Jac_down_trunc_132= 0 ; %16
Jac_down_trunc_232= 0 ; %17
Jac_down_trunc_332= 0 ; %18
Jac_down_trunc_113= 0 ; %19
Jac_down_trunc_213= 0 ; %20
Jac_down_trunc_313= 0 ; %21
Jac_down_trunc_123= 0 ; %22
Jac_down_trunc_223= 0 ; %23
Jac_down_trunc_323= 0 ; %24
Jac_down_trunc_133= 0 ; %25
Jac_down_trunc_233= 0 ; %26
Jac_down_trunc_333= 0 ; %27

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

Flux_x = @(q,p)[ q(2);
    q(2)^2/q(1) + p;
 (q(2)*q(3))/q(1);
 (q(2)*q(4))/q(1);
 (q(5)*q(2))/q(1);
 (q(6)*q(2))/q(1);
 (q(7)*q(2))/q(1);];
 
Flux_y = @(q,p)[ q(3);
 (q(2)*q(3))/q(1);
    q(3)^2/q(1) + p ;
 (q(3)*q(4))/q(1);
 (q(5)*q(3))/q(1);
 (q(6)*q(3))/q(1);
 (q(7)*q(3))/q(1);];
 
Flux_z = @(q,p)[ q(4);
 (q(2)*q(4))/q(1);
 (q(3)*q(4))/q(1);
    q(4)^2/q(1) + p ;
 (q(5)*q(4))/q(1);
 (q(6)*q(4))/q(1);
 (q(7)*q(4))/q(1);];


DFlux_dq_x = @(q,dpdn) [[               0,        1,    0,    0,    0,    0,    0];
[ dpdn - q(2)^2/q(1)^2, (2*q(2))/q(1),    0,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2,     q(3)/q(1), q(2)/q(1),    0,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2,     q(4)/q(1),    0, q(2)/q(1),    0,    0,    0];
[    -(q(5)*q(2))/q(1)^2,     q(5)/q(1),    0,    0, q(2)/q(1),    0,    0];
[    -(q(6)*q(2))/q(1)^2,     q(6)/q(1),    0,    0,    0, q(2)/q(1),    0];
[    -(q(7)*q(2))/q(1)^2,     q(7)/q(1),    0,    0,    0,    0, q(2)/q(1)]];


DFlux_dq_y = @(q,dpdn) [[               0,    0,        1,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2, q(3)/q(1),     q(2)/q(1),    0,    0,    0,    0];
[ dpdn - q(3)^2/q(1)^2,    0, (2*q(3))/q(1),    0,    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0,     q(4)/q(1), q(3)/q(1),    0,    0,    0];
[    -(q(5)*q(3))/q(1)^2,    0,     q(5)/q(1),    0, q(3)/q(1),    0,    0];
[    -(q(6)*q(3))/q(1)^2,    0,     q(6)/q(1),    0,    0, q(3)/q(1),    0];
[    -(q(7)*q(3))/q(1)^2,    0,     q(7)/q(1),    0,    0,    0, q(3)/q(1)]];

DFlux_dq_z = @(q,dpdn) [[               0,    0,    0,        1,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2, q(4)/q(1),    0,     q(2)/q(1),    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0, q(4)/q(1),     q(3)/q(1),    0,    0,    0];
[ dpdn - q(4)^2/q(1)^2,    0,    0, (2*q(4))/q(1),    0,    0,    0];
[    -(q(5)*q(4))/q(1)^2,    0,    0,     q(5)/q(1), q(4)/q(1),    0,    0];
[    -(q(6)*q(4))/q(1)^2,    0,    0,     q(6)/q(1),    0, q(4)/q(1),    0];
[    -(q(7)*q(4))/q(1)^2,    0,    0,     q(7)/q(1),    0,    0, q(4)/q(1)];];


speed_x = @(q,dpdn) abs(q(2)/q(1)) + abs(sqrt(dpdn)) ; 
speed_y = @(q,dpdn) abs(q(3)/q(1)) + abs(sqrt(dpdn)) ; 
speed_z = @(q,dpdn) abs(q(4)/q(1)) + abs(sqrt(dpdn)) ; 


for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = uadlocs(k,1);
    v1quad = Dlist(k,2);
    v2quad = Dlist(k,3);
    
    weight = Dquadwgts(k);    
    wgtnuv1 = weight*nuv1;
    wgtnuv2 = weight*nuv2;
    wgtnuv3 = weight*nuv3;
    
    psie = vectpsi_east(:,:,tquad,v1quad,v2quad);
    psiw = vectpsi_west(:,:,tquad,v1quad,v2quad);
    psin = vectpsi_nort(:,:,tquad,v1quad,v2quad);
    psis = vectpsi_sout(:,:,tquad,v1quad,v2quad);
    psiu = vectpsi_uppr(:,:,tquad,v1quad,v2quad);
    psid = vectpsi_down(:,:,tquad,v1quad,v2quad);
    
    lilphi_e = vectphi_east(1,1:data.Ls,v1quad,v2quad);
    lilphi_w = vectphi_west(1,1:data.Ls,v1quad,v2quad);
    lilphi_n = vectphi_nort(1,1:data.Ls,v1quad,v2quad);
    lilphi_s = vectphi_sout(1,1:data.Ls,v1quad,v2quad);
    lilphi_u = vectphi_uppr(1,1:data.Ls,v1quad,v2quad);
    lilphi_d = vectphi_down(1,1:data.Ls,v1quad,v2quad);
    
    psieT = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    psiwT = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    psinT = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    psisT = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    psiuT = vectpsi_uppr_Trans(:,:,tquad,v1quad,v2quad);
    psidT = vectpsi_down_Trans(:,:,tquad,v1quad,v2quad);
    
    % EAST AND WEST TERMS
    w01_e = psie*q01;
    w01_w = psiw*q01;
    
    w02_e = psie*q02;
    w02_w = psiw*q02;
    
    w03_e = psie*q03;
    w03_w = psiw*q03;    
    
    w04_e = psie*q04;
    w04_w = psiw*q04;    

    w05_e = psie*q05;
    w05_w = psiw*q05;    

    w06_e = psie*q06;
    w06_w = psiw*q06;    

    w07_e = psie*q07;
    w07_w = psiw*q07;    

    w08_e = psie*q08;
    w08_w = psiw*q08;    

    w09_e = psie*q09;
    w09_w = psiw*q09;
    
    w10_e = psie*q10;
    w10_w = psiw*q10;
    
    w11_e = psie*q11;
    w11_w = psiw*q11;

    w12_e = psie*q12;
    w12_w = psiw*q12;
   
    w13_e = psie*q13;
    w13_w = psiw*q13;    
    
    w14_e = psie*q14;
    w14_w = psiw*q14;    

    w15_e = psie*q15;
    w15_w = psiw*q15;    

    w16_e = psie*q16;
    w16_w = psiw*q16;    

    w17_e = psie*q17;
    w17_w = psiw*q17;    

    w18_e = psie*q18;
    w18_w = psiw*q18;    

    w19_e = psie*q19;
    w19_w = psiw*q19;
    
    w20_e = psie*q20;
    w20_w = psiw*q20;
    
    w21_e = psie*q21;
    w21_w = psiw*q21;

    w22_e = psie*q22;
    w22_w = psiw*q22;
   
    w23_e = psie*q23;
    w23_w = psiw*q23;    
    
    w24_e = psie*q24;
    w24_w = psiw*q24;    

    w25_e = psie*q25;
    w25_w = psiw*q25;    

    w26_e = psie*q26;
    w26_w = psiw*q26;    

    w27_e = psie*q27;
    w27_w = psiw*q27;      

% NORT SOUT TERMS
    w01_n = psin*q01;
    w01_s = psis*q01;
    
    w02_n = psin*q02;
    w02_s = psis*q02;
    
    w03_n = psin*q03;
    w03_s = psis*q03;    
    
    w04_n = psin*q04;
    w04_s = psis*q04;    

    w05_n = psin*q05;
    w05_s = psis*q05;    

    w06_n = psin*q06;
    w06_s = psis*q06;    

    w07_n = psin*q07;
    w07_s = psis*q07;    

    w08_n = psin*q08;
    w08_s = psis*q08;    

    w09_n = psin*q09;
    w09_s = psis*q09;
    
    w10_n = psin*q10;
    w10_s = psis*q10;
    
    w11_n = psin*q11;
    w11_s = psis*q11;

    w12_n = psin*q12;
    w12_s = psis*q12;
   
    w13_n = psin*q13;
    w13_s = psis*q13;    
    
    w14_n = psin*q14;
    w14_s = psis*q14;    

    w15_n = psin*q15;
    w15_s = psis*q15;    

    w16_n = psin*q16;
    w16_s = psis*q16;    

    w17_n = psin*q17;
    w17_s = psis*q17;    

    w18_n = psin*q18;
    w18_s = psis*q18;    

    w19_n = psin*q19;
    w19_s = psis*q19;
    
    w20_n = psin*q20;
    w20_s = psis*q20;
    
    w21_n = psin*q21;
    w21_s = psis*q21;

    w22_n = psin*q22;
    w22_s = psis*q22;
   
    w23_n = psin*q23;
    w23_s = psis*q23;    
    
    w24_n = psin*q24;
    w24_s = psis*q24;    

    w25_n = psin*q25;
    w25_s = psis*q25;    

    w26_n = psin*q26;
    w26_s = psis*q26;    

    w27_n = psin*q27;
    w27_s = psis*q27;  
    
   

% UPPR DOWN TERMS
    w01_u = psiu*q01;
    w01_d = psid*q01;
    
    w02_u = psiu*q02;
    w02_d = psid*q02;
    
    w03_u = psiu*q03;
    w03_d = psid*q03;    
    
    w04_u = psiu*q04;
    w04_d = psid*q04;    

    w05_u = psiu*q05;
    w05_d = psid*q05;    

    w06_u = psiu*q06;
    w06_d = psid*q06;    

    w07_u = psiu*q07;
    w07_d = psid*q07;    

    w08_u = psiu*q08;
    w08_d = psid*q08;    

    w09_u = psiu*q09;
    w09_d = psid*q09;
    
    w10_u = psiu*q10;
    w10_d = psid*q10;
    
    w11_u = psiu*q11;
    w11_d = psid*q11;

    w12_u = psiu*q12;
    w12_d = psid*q12;
   
    w13_u = psiu*q13;
    w13_d = psid*q13;    
    
    w14_u = psiu*q14;
    w14_d = psid*q14;    

    w15_u = psiu*q15;
    w15_d = psid*q15;    

    w16_u = psiu*q16;
    w16_d = psid*q16;    

    w17_u = psiu*q17;
    w17_d = psid*q17;    

    w18_u = psiu*q18;
    w18_d = psid*q18;    

    w19_u = psiu*q19;
    w19_d = psid*q19;
    
    w20_u = psiu*q20;
    w20_d = psid*q20;
    
    w21_u = psiu*q21;
    w21_d = psid*q21;

    w22_u = psiu*q22;
    w22_d = psid*q22;
   
    w23_u = psiu*q23;
    w23_d = psid*q23;    
    
    w24_u = psiu*q24;
    w24_d = psid*q24;    

    w25_u = psiu*q25;
    w25_d = psid*q25;    

    w26_u = psiu*q26;
    w26_d = psid*q26;    

    w27_u = psiu*q27;
    w27_d = psid*q27;  
    
    p01_e= lilphi_e*P01;
    p02_e= lilphi_e*P02;
    p03_e= lilphi_e*P03;
    p04_e= lilphi_e*P04;
    p05_e= lilphi_e*P05;
    p06_e= lilphi_e*P06;
    p07_e= lilphi_e*P07;
    p08_e= lilphi_e*P08;
    p09_e= lilphi_e*P09;
    p10_e= lilphi_e*P10;
    p11_e= lilphi_e*P11;
    p12_e= lilphi_e*P12;
    p13_e= lilphi_e*P13;
    p14_e= lilphi_e*P14;
    p15_e= lilphi_e*P15;
    p16_e= lilphi_e*P16;
    p17_e= lilphi_e*P17;
    p18_e= lilphi_e*P18;
    p19_e= lilphi_e*P19;
    p20_e= lilphi_e*P20;
    p21_e= lilphi_e*P21;
    p22_e= lilphi_e*P22;
    p23_e= lilphi_e*P23;
    p24_e= lilphi_e*P24;
    p25_e= lilphi_e*P25;
    p26_e= lilphi_e*P26;
    p27_e= lilphi_e*P27;
    p01_w= lilphi_w*P01;
    p02_w= lilphi_w*P02;
    p03_w= lilphi_w*P03;
    p04_w= lilphi_w*P04;
    p05_w= lilphi_w*P05;
    p06_w= lilphi_w*P06;
    p07_w= lilphi_w*P07;
    p08_w= lilphi_w*P08;
    p09_w= lilphi_w*P09;
    p10_w= lilphi_w*P10;
    p11_w= lilphi_w*P11;
    p12_w= lilphi_w*P12;
    p13_w= lilphi_w*P13;
    p14_w= lilphi_w*P14;
    p15_w= lilphi_w*P15;
    p16_w= lilphi_w*P16;
    p17_w= lilphi_w*P17;
    p18_w= lilphi_w*P18;
    p19_w= lilphi_w*P19;
    p20_w= lilphi_w*P20;
    p21_w= lilphi_w*P21;
    p22_w= lilphi_w*P22;
    p23_w= lilphi_w*P23;
    p24_w= lilphi_w*P24;
    p25_w= lilphi_w*P25;
    p26_w= lilphi_w*P26;
    p27_w= lilphi_w*P27;
    p01_n= lilphi_n*P01;
    p02_n= lilphi_n*P02;
    p03_n= lilphi_n*P03;
    p04_n= lilphi_n*P04;
    p05_n= lilphi_n*P05;
    p06_n= lilphi_n*P06;
    p07_n= lilphi_n*P07;
    p08_n= lilphi_n*P08;
    p09_n= lilphi_n*P09;
    p10_n= lilphi_n*P10;
    p11_n= lilphi_n*P11;
    p12_n= lilphi_n*P12;
    p13_n= lilphi_n*P13;
    p14_n= lilphi_n*P14;
    p15_n= lilphi_n*P15;
    p16_n= lilphi_n*P16;
    p17_n= lilphi_n*P17;
    p18_n= lilphi_n*P18;
    p19_n= lilphi_n*P19;
    p20_n= lilphi_n*P20;
    p21_n= lilphi_n*P21;
    p22_n= lilphi_n*P22;
    p23_n= lilphi_n*P23;
    p24_n= lilphi_n*P24;
    p25_n= lilphi_n*P25;
    p26_n= lilphi_n*P26;
    p27_n= lilphi_n*P27;
    p01_s= lilphi_s*P01;
    p02_s= lilphi_s*P02;
    p03_s= lilphi_s*P03;
    p04_s= lilphi_s*P04;
    p05_s= lilphi_s*P05;
    p06_s= lilphi_s*P06;
    p07_s= lilphi_s*P07;
    p08_s= lilphi_s*P08;
    p09_s= lilphi_s*P09;
    p10_s= lilphi_s*P10;
    p11_s= lilphi_s*P11;
    p12_s= lilphi_s*P12;
    p13_s= lilphi_s*P13;
    p14_s= lilphi_s*P14;
    p15_s= lilphi_s*P15;
    p16_s= lilphi_s*P16;
    p17_s= lilphi_s*P17;
    p18_s= lilphi_s*P18;
    p19_s= lilphi_s*P19;
    p20_s= lilphi_s*P20;
    p21_s= lilphi_s*P21;
    p22_s= lilphi_s*P22;
    p23_s= lilphi_s*P23;
    p24_s= lilphi_s*P24;
    p25_s= lilphi_s*P25;
    p26_s= lilphi_s*P26;
    p27_s= lilphi_s*P27;
    p01_u= lilphi_u*P01;
    p02_u= lilphi_u*P02;
    p03_u= lilphi_u*P03;
    p04_u= lilphi_u*P04;
    p05_u= lilphi_u*P05;
    p06_u= lilphi_u*P06;
    p07_u= lilphi_u*P07;
    p08_u= lilphi_u*P08;
    p09_u= lilphi_u*P09;
    p10_u= lilphi_u*P10;
    p11_u= lilphi_u*P11;
    p12_u= lilphi_u*P12;
    p13_u= lilphi_u*P13;
    p14_u= lilphi_u*P14;
    p15_u= lilphi_u*P15;
    p16_u= lilphi_u*P16;
    p17_u= lilphi_u*P17;
    p18_u= lilphi_u*P18;
    p19_u= lilphi_u*P19;
    p20_u= lilphi_u*P20;
    p21_u= lilphi_u*P21;
    p22_u= lilphi_u*P22;
    p23_u= lilphi_u*P23;
    p24_u= lilphi_u*P24;
    p25_u= lilphi_u*P25;
    p26_u= lilphi_u*P26;
    p27_u= lilphi_u*P27;
    p01_d= lilphi_d*P01;
    p02_d= lilphi_d*P02;
    p03_d= lilphi_d*P03;
    p04_d= lilphi_d*P04;
    p05_d= lilphi_d*P05;
    p06_d= lilphi_d*P06;
    p07_d= lilphi_d*P07;
    p08_d= lilphi_d*P08;
    p09_d= lilphi_d*P09;
    p10_d= lilphi_d*P10;
    p11_d= lilphi_d*P11;
    p12_d= lilphi_d*P12;
    p13_d= lilphi_d*P13;
    p14_d= lilphi_d*P14;
    p15_d= lilphi_d*P15;
    p16_d= lilphi_d*P16;
    p17_d= lilphi_d*P17;
    p18_d= lilphi_d*P18;
    p19_d= lilphi_d*P19;
    p20_d= lilphi_d*P20;
    p21_d= lilphi_d*P21;
    p22_d= lilphi_d*P22;
    p23_d= lilphi_d*P23;
    p24_d= lilphi_d*P24;
    p25_d= lilphi_d*P25;
    p26_d= lilphi_d*P26;
    p27_d= lilphi_d*P27;    
    
    dpdn01_e= lilphi_e*DPDN01;
    dpdn02_e= lilphi_e*DPDN02;
    dpdn03_e= lilphi_e*DPDN03;
    dpdn04_e= lilphi_e*DPDN04;
    dpdn05_e= lilphi_e*DPDN05;
    dpdn06_e= lilphi_e*DPDN06;
    dpdn07_e= lilphi_e*DPDN07;
    dpdn08_e= lilphi_e*DPDN08;
    dpdn09_e= lilphi_e*DPDN09;
    dpdn10_e= lilphi_e*DPDN10;
    dpdn11_e= lilphi_e*DPDN11;
    dpdn12_e= lilphi_e*DPDN12;
    dpdn13_e= lilphi_e*DPDN13;
    dpdn14_e= lilphi_e*DPDN14;
    dpdn15_e= lilphi_e*DPDN15;
    dpdn16_e= lilphi_e*DPDN16;
    dpdn17_e= lilphi_e*DPDN17;
    dpdn18_e= lilphi_e*DPDN18;
    dpdn19_e= lilphi_e*DPDN19;
    dpdn20_e= lilphi_e*DPDN20;
    dpdn21_e= lilphi_e*DPDN21;
    dpdn22_e= lilphi_e*DPDN22;
    dpdn23_e= lilphi_e*DPDN23;
    dpdn24_e= lilphi_e*DPDN24;
    dpdn25_e= lilphi_e*DPDN25;
    dpdn26_e= lilphi_e*DPDN26;
    dpdn27_e= lilphi_e*DPDN27;
    dpdn01_w= lilphi_w*DPDN01;
    dpdn02_w= lilphi_w*DPDN02;
    dpdn03_w= lilphi_w*DPDN03;
    dpdn04_w= lilphi_w*DPDN04;
    dpdn05_w= lilphi_w*DPDN05;
    dpdn06_w= lilphi_w*DPDN06;
    dpdn07_w= lilphi_w*DPDN07;
    dpdn08_w= lilphi_w*DPDN08;
    dpdn09_w= lilphi_w*DPDN09;
    dpdn10_w= lilphi_w*DPDN10;
    dpdn11_w= lilphi_w*DPDN11;
    dpdn12_w= lilphi_w*DPDN12;
    dpdn13_w= lilphi_w*DPDN13;
    dpdn14_w= lilphi_w*DPDN14;
    dpdn15_w= lilphi_w*DPDN15;
    dpdn16_w= lilphi_w*DPDN16;
    dpdn17_w= lilphi_w*DPDN17;
    dpdn18_w= lilphi_w*DPDN18;
    dpdn19_w= lilphi_w*DPDN19;
    dpdn20_w= lilphi_w*DPDN20;
    dpdn21_w= lilphi_w*DPDN21;
    dpdn22_w= lilphi_w*DPDN22;
    dpdn23_w= lilphi_w*DPDN23;
    dpdn24_w= lilphi_w*DPDN24;
    dpdn25_w= lilphi_w*DPDN25;
    dpdn26_w= lilphi_w*DPDN26;
    dpdn27_w= lilphi_w*DPDN27;
    dpdn01_n= lilphi_n*DPDN01;
    dpdn02_n= lilphi_n*DPDN02;
    dpdn03_n= lilphi_n*DPDN03;
    dpdn04_n= lilphi_n*DPDN04;
    dpdn05_n= lilphi_n*DPDN05;
    dpdn06_n= lilphi_n*DPDN06;
    dpdn07_n= lilphi_n*DPDN07;
    dpdn08_n= lilphi_n*DPDN08;
    dpdn09_n= lilphi_n*DPDN09;
    dpdn10_n= lilphi_n*DPDN10;
    dpdn11_n= lilphi_n*DPDN11;
    dpdn12_n= lilphi_n*DPDN12;
    dpdn13_n= lilphi_n*DPDN13;
    dpdn14_n= lilphi_n*DPDN14;
    dpdn15_n= lilphi_n*DPDN15;
    dpdn16_n= lilphi_n*DPDN16;
    dpdn17_n= lilphi_n*DPDN17;
    dpdn18_n= lilphi_n*DPDN18;
    dpdn19_n= lilphi_n*DPDN19;
    dpdn20_n= lilphi_n*DPDN20;
    dpdn21_n= lilphi_n*DPDN21;
    dpdn22_n= lilphi_n*DPDN22;
    dpdn23_n= lilphi_n*DPDN23;
    dpdn24_n= lilphi_n*DPDN24;
    dpdn25_n= lilphi_n*DPDN25;
    dpdn26_n= lilphi_n*DPDN26;
    dpdn27_n= lilphi_n*DPDN27;
    dpdn01_s= lilphi_s*DPDN01;
    dpdn02_s= lilphi_s*DPDN02;
    dpdn03_s= lilphi_s*DPDN03;
    dpdn04_s= lilphi_s*DPDN04;
    dpdn05_s= lilphi_s*DPDN05;
    dpdn06_s= lilphi_s*DPDN06;
    dpdn07_s= lilphi_s*DPDN07;
    dpdn08_s= lilphi_s*DPDN08;
    dpdn09_s= lilphi_s*DPDN09;
    dpdn10_s= lilphi_s*DPDN10;
    dpdn11_s= lilphi_s*DPDN11;
    dpdn12_s= lilphi_s*DPDN12;
    dpdn13_s= lilphi_s*DPDN13;
    dpdn14_s= lilphi_s*DPDN14;
    dpdn15_s= lilphi_s*DPDN15;
    dpdn16_s= lilphi_s*DPDN16;
    dpdn17_s= lilphi_s*DPDN17;
    dpdn18_s= lilphi_s*DPDN18;
    dpdn19_s= lilphi_s*DPDN19;
    dpdn20_s= lilphi_s*DPDN20;
    dpdn21_s= lilphi_s*DPDN21;
    dpdn22_s= lilphi_s*DPDN22;
    dpdn23_s= lilphi_s*DPDN23;
    dpdn24_s= lilphi_s*DPDN24;
    dpdn25_s= lilphi_s*DPDN25;
    dpdn26_s= lilphi_s*DPDN26;
    dpdn27_s= lilphi_s*DPDN27;
    dpdn01_u= lilphi_u*DPDN01;
    dpdn02_u= lilphi_u*DPDN02;
    dpdn03_u= lilphi_u*DPDN03;
    dpdn04_u= lilphi_u*DPDN04;
    dpdn05_u= lilphi_u*DPDN05;
    dpdn06_u= lilphi_u*DPDN06;
    dpdn07_u= lilphi_u*DPDN07;
    dpdn08_u= lilphi_u*DPDN08;
    dpdn09_u= lilphi_u*DPDN09;
    dpdn10_u= lilphi_u*DPDN10;
    dpdn11_u= lilphi_u*DPDN11;
    dpdn12_u= lilphi_u*DPDN12;
    dpdn13_u= lilphi_u*DPDN13;
    dpdn14_u= lilphi_u*DPDN14;
    dpdn15_u= lilphi_u*DPDN15;
    dpdn16_u= lilphi_u*DPDN16;
    dpdn17_u= lilphi_u*DPDN17;
    dpdn18_u= lilphi_u*DPDN18;
    dpdn19_u= lilphi_u*DPDN19;
    dpdn20_u= lilphi_u*DPDN20;
    dpdn21_u= lilphi_u*DPDN21;
    dpdn22_u= lilphi_u*DPDN22;
    dpdn23_u= lilphi_u*DPDN23;
    dpdn24_u= lilphi_u*DPDN24;
    dpdn25_u= lilphi_u*DPDN25;
    dpdn26_u= lilphi_u*DPDN26;
    dpdn27_u= lilphi_u*DPDN27;
    dpdn01_d= lilphi_d*DPDN01;
    dpdn02_d= lilphi_d*DPDN02;
    dpdn03_d= lilphi_d*DPDN03;
    dpdn04_d= lilphi_d*DPDN04;
    dpdn05_d= lilphi_d*DPDN05;
    dpdn06_d= lilphi_d*DPDN06;
    dpdn07_d= lilphi_d*DPDN07;
    dpdn08_d= lilphi_d*DPDN08;
    dpdn09_d= lilphi_d*DPDN09;
    dpdn10_d= lilphi_d*DPDN10;
    dpdn11_d= lilphi_d*DPDN11;
    dpdn12_d= lilphi_d*DPDN12;
    dpdn13_d= lilphi_d*DPDN13;
    dpdn14_d= lilphi_d*DPDN14;
    dpdn15_d= lilphi_d*DPDN15;
    dpdn16_d= lilphi_d*DPDN16;
    dpdn17_d= lilphi_d*DPDN17;
    dpdn18_d= lilphi_d*DPDN18;
    dpdn19_d= lilphi_d*DPDN19;
    dpdn20_d= lilphi_d*DPDN20;
    dpdn21_d= lilphi_d*DPDN21;
    dpdn22_d= lilphi_d*DPDN22;
    dpdn23_d= lilphi_d*DPDN23;
    dpdn24_d= lilphi_d*DPDN24;
    dpdn25_d= lilphi_d*DPDN25;
    dpdn26_d= lilphi_d*DPDN26;
    dpdn27_d= lilphi_d*DPDN27;
      
    F01_e= Flux_x(w01_e,p01_e);
    F02_e= Flux_x(w02_e,p02_e);
    F03_e= Flux_x(w03_e,p03_e);
    F04_e= Flux_x(w04_e,p04_e);
    F05_e= Flux_x(w05_e,p05_e);
    F06_e= Flux_x(w06_e,p06_e);
    F07_e= Flux_x(w07_e,p07_e);
    F08_e= Flux_x(w08_e,p08_e);
    F09_e= Flux_x(w09_e,p09_e);
    F10_e= Flux_x(w10_e,p10_e);
    F11_e= Flux_x(w11_e,p11_e);
    F12_e= Flux_x(w12_e,p12_e);
    F13_e= Flux_x(w13_e,p13_e);
    F14_e= Flux_x(w14_e,p14_e);
    F15_e= Flux_x(w15_e,p15_e);
    F16_e= Flux_x(w16_e,p16_e);
    F17_e= Flux_x(w17_e,p17_e);
    F18_e= Flux_x(w18_e,p18_e);
    F19_e= Flux_x(w19_e,p19_e);
    F20_e= Flux_x(w20_e,p20_e);
    F21_e= Flux_x(w21_e,p21_e);
    F22_e= Flux_x(w22_e,p22_e);
    F23_e= Flux_x(w23_e,p23_e);
    F24_e= Flux_x(w24_e,p24_e);
    F25_e= Flux_x(w25_e,p25_e);
    F26_e= Flux_x(w26_e,p26_e);
    F27_e= Flux_x(w27_e,p27_e);
    F01_w= Flux_x(w01_w,p01_w);
    F02_w= Flux_x(w02_w,p02_w);
    F03_w= Flux_x(w03_w,p03_w);
    F04_w= Flux_x(w04_w,p04_w);
    F05_w= Flux_x(w05_w,p05_w);
    F06_w= Flux_x(w06_w,p06_w);
    F07_w= Flux_x(w07_w,p07_w);
    F08_w= Flux_x(w08_w,p08_w);
    F09_w= Flux_x(w09_w,p09_w);
    F10_w= Flux_x(w10_w,p10_w);
    F11_w= Flux_x(w11_w,p11_w);
    F12_w= Flux_x(w12_w,p12_w);
    F13_w= Flux_x(w13_w,p13_w);
    F14_w= Flux_x(w14_w,p14_w);
    F15_w= Flux_x(w15_w,p15_w);
    F16_w= Flux_x(w16_w,p16_w);
    F17_w= Flux_x(w17_w,p17_w);
    F18_w= Flux_x(w18_w,p18_w);
    F19_w= Flux_x(w19_w,p19_w);
    F20_w= Flux_x(w20_w,p20_w);
    F21_w= Flux_x(w21_w,p21_w);
    F22_w= Flux_x(w22_w,p22_w);
    F23_w= Flux_x(w23_w,p23_w);
    F24_w= Flux_x(w24_w,p24_w);
    F25_w= Flux_x(w25_w,p25_w);
    F26_w= Flux_x(w26_w,p26_w);
    F27_w= Flux_x(w27_w,p27_w);
    F01_n= Flux_y(w01_n,p01_n);
    F02_n= Flux_y(w02_n,p02_n);
    F03_n= Flux_y(w03_n,p03_n);
    F04_n= Flux_y(w04_n,p04_n);
    F05_n= Flux_y(w05_n,p05_n);
    F06_n= Flux_y(w06_n,p06_n);
    F07_n= Flux_y(w07_n,p07_n);
    F08_n= Flux_y(w08_n,p08_n);
    F09_n= Flux_y(w09_n,p09_n);
    F10_n= Flux_y(w10_n,p10_n);
    F11_n= Flux_y(w11_n,p11_n);
    F12_n= Flux_y(w12_n,p12_n);
    F13_n= Flux_y(w13_n,p13_n);
    F14_n= Flux_y(w14_n,p14_n);
    F15_n= Flux_y(w15_n,p15_n);
    F16_n= Flux_y(w16_n,p16_n);
    F17_n= Flux_y(w17_n,p17_n);
    F18_n= Flux_y(w18_n,p18_n);
    F19_n= Flux_y(w19_n,p19_n);
    F20_n= Flux_y(w20_n,p20_n);
    F21_n= Flux_y(w21_n,p21_n);
    F22_n= Flux_y(w22_n,p22_n);
    F23_n= Flux_y(w23_n,p23_n);
    F24_n= Flux_y(w24_n,p24_n);
    F25_n= Flux_y(w25_n,p25_n);
    F26_n= Flux_y(w26_n,p26_n);
    F27_n= Flux_y(w27_n,p27_n);
    F01_s= Flux_y(w01_s,p01_s);
    F02_s= Flux_y(w02_s,p02_s);
    F03_s= Flux_y(w03_s,p03_s);
    F04_s= Flux_y(w04_s,p04_s);
    F05_s= Flux_y(w05_s,p05_s);
    F06_s= Flux_y(w06_s,p06_s);
    F07_s= Flux_y(w07_s,p07_s);
    F08_s= Flux_y(w08_s,p08_s);
    F09_s= Flux_y(w09_s,p09_s);
    F10_s= Flux_y(w10_s,p10_s);
    F11_s= Flux_y(w11_s,p11_s);
    F12_s= Flux_y(w12_s,p12_s);
    F13_s= Flux_y(w13_s,p13_s);
    F14_s= Flux_y(w14_s,p14_s);
    F15_s= Flux_y(w15_s,p15_s);
    F16_s= Flux_y(w16_s,p16_s);
    F17_s= Flux_y(w17_s,p17_s);
    F18_s= Flux_y(w18_s,p18_s);
    F19_s= Flux_y(w19_s,p19_s);
    F20_s= Flux_y(w20_s,p20_s);
    F21_s= Flux_y(w21_s,p21_s);
    F22_s= Flux_y(w22_s,p22_s);
    F23_s= Flux_y(w23_s,p23_s);
    F24_s= Flux_y(w24_s,p24_s);
    F25_s= Flux_y(w25_s,p25_s);
    F26_s= Flux_y(w26_s,p26_s);
    F27_s= Flux_y(w27_s,p27_s);
    F01_u= Flux_z(w01_u,p01_u);
    F02_u= Flux_z(w02_u,p02_u);
    F03_u= Flux_z(w03_u,p03_u);
    F04_u= Flux_z(w04_u,p04_u);
    F05_u= Flux_z(w05_u,p05_u);
    F06_u= Flux_z(w06_u,p06_u);
    F07_u= Flux_z(w07_u,p07_u);
    F08_u= Flux_z(w08_u,p08_u);
    F09_u= Flux_z(w09_u,p09_u);
    F10_u= Flux_z(w10_u,p10_u);
    F11_u= Flux_z(w11_u,p11_u);
    F12_u= Flux_z(w12_u,p12_u);
    F13_u= Flux_z(w13_u,p13_u);
    F14_u= Flux_z(w14_u,p14_u);
    F15_u= Flux_z(w15_u,p15_u);
    F16_u= Flux_z(w16_u,p16_u);
    F17_u= Flux_z(w17_u,p17_u);
    F18_u= Flux_z(w18_u,p18_u);
    F19_u= Flux_z(w19_u,p19_u);
    F20_u= Flux_z(w20_u,p20_u);
    F21_u= Flux_z(w21_u,p21_u);
    F22_u= Flux_z(w22_u,p22_u);
    F23_u= Flux_z(w23_u,p23_u);
    F24_u= Flux_z(w24_u,p24_u);
    F25_u= Flux_z(w25_u,p25_u);
    F26_u= Flux_z(w26_u,p26_u);
    F27_u= Flux_z(w27_u,p27_u);
    F01_d= Flux_z(w01_d,p01_d);
    F02_d= Flux_z(w02_d,p02_d);
    F03_d= Flux_z(w03_d,p03_d);
    F04_d= Flux_z(w04_d,p04_d);
    F05_d= Flux_z(w05_d,p05_d);
    F06_d= Flux_z(w06_d,p06_d);
    F07_d= Flux_z(w07_d,p07_d);
    F08_d= Flux_z(w08_d,p08_d);
    F09_d= Flux_z(w09_d,p09_d);
    F10_d= Flux_z(w10_d,p10_d);
    F11_d= Flux_z(w11_d,p11_d);
    F12_d= Flux_z(w12_d,p12_d);
    F13_d= Flux_z(w13_d,p13_d);
    F14_d= Flux_z(w14_d,p14_d);
    F15_d= Flux_z(w15_d,p15_d);
    F16_d= Flux_z(w16_d,p16_d);
    F17_d= Flux_z(w17_d,p17_d);
    F18_d= Flux_z(w18_d,p18_d);
    F19_d= Flux_z(w19_d,p19_d);
    F20_d= Flux_z(w20_d,p20_d);
    F21_d= Flux_z(w21_d,p21_d);
    F22_d= Flux_z(w22_d,p22_d);
    F23_d= Flux_z(w23_d,p23_d);
    F24_d= Flux_z(w24_d,p24_d);
    F25_d= Flux_z(w25_d,p25_d);
    F26_d= Flux_z(w26_d,p26_d);
    F27_d= Flux_z(w27_d,p27_d);

    DF01_e= DFlux_dq_x(w01_e,dpdn01_e);
    DF02_e= DFlux_dq_x(w02_e,dpdn02_e);
    DF03_e= DFlux_dq_x(w03_e,dpdn03_e);
    DF04_e= DFlux_dq_x(w04_e,dpdn04_e);
    DF05_e= DFlux_dq_x(w05_e,dpdn05_e);
    DF06_e= DFlux_dq_x(w06_e,dpdn06_e);
    DF07_e= DFlux_dq_x(w07_e,dpdn07_e);
    DF08_e= DFlux_dq_x(w08_e,dpdn08_e);
    DF09_e= DFlux_dq_x(w09_e,dpdn09_e);
    DF10_e= DFlux_dq_x(w10_e,dpdn10_e);
    DF11_e= DFlux_dq_x(w11_e,dpdn11_e);
    DF12_e= DFlux_dq_x(w12_e,dpdn12_e);
    DF13_e= DFlux_dq_x(w13_e,dpdn13_e);
    DF14_e= DFlux_dq_x(w14_e,dpdn14_e);
    DF15_e= DFlux_dq_x(w15_e,dpdn15_e);
    DF16_e= DFlux_dq_x(w16_e,dpdn16_e);
    DF17_e= DFlux_dq_x(w17_e,dpdn17_e);
    DF18_e= DFlux_dq_x(w18_e,dpdn18_e);
    DF19_e= DFlux_dq_x(w19_e,dpdn19_e);
    DF20_e= DFlux_dq_x(w20_e,dpdn20_e);
    DF21_e= DFlux_dq_x(w21_e,dpdn21_e);
    DF22_e= DFlux_dq_x(w22_e,dpdn22_e);
    DF23_e= DFlux_dq_x(w23_e,dpdn23_e);
    DF24_e= DFlux_dq_x(w24_e,dpdn24_e);
    DF25_e= DFlux_dq_x(w25_e,dpdn25_e);
    DF26_e= DFlux_dq_x(w26_e,dpdn26_e);
    DF27_e= DFlux_dq_x(w27_e,dpdn27_e);
    DF01_w= DFlux_dq_x(w01_w,dpdn01_w);
    DF02_w= DFlux_dq_x(w02_w,dpdn02_w);
    DF03_w= DFlux_dq_x(w03_w,dpdn03_w);
    DF04_w= DFlux_dq_x(w04_w,dpdn04_w);
    DF05_w= DFlux_dq_x(w05_w,dpdn05_w);
    DF06_w= DFlux_dq_x(w06_w,dpdn06_w);
    DF07_w= DFlux_dq_x(w07_w,dpdn07_w);
    DF08_w= DFlux_dq_x(w08_w,dpdn08_w);
    DF09_w= DFlux_dq_x(w09_w,dpdn09_w);
    DF10_w= DFlux_dq_x(w10_w,dpdn10_w);
    DF11_w= DFlux_dq_x(w11_w,dpdn11_w);
    DF12_w= DFlux_dq_x(w12_w,dpdn12_w);
    DF13_w= DFlux_dq_x(w13_w,dpdn13_w);
    DF14_w= DFlux_dq_x(w14_w,dpdn14_w);
    DF15_w= DFlux_dq_x(w15_w,dpdn15_w);
    DF16_w= DFlux_dq_x(w16_w,dpdn16_w);
    DF17_w= DFlux_dq_x(w17_w,dpdn17_w);
    DF18_w= DFlux_dq_x(w18_w,dpdn18_w);
    DF19_w= DFlux_dq_x(w19_w,dpdn19_w);
    DF20_w= DFlux_dq_x(w20_w,dpdn20_w);
    DF21_w= DFlux_dq_x(w21_w,dpdn21_w);
    DF22_w= DFlux_dq_x(w22_w,dpdn22_w);
    DF23_w= DFlux_dq_x(w23_w,dpdn23_w);
    DF24_w= DFlux_dq_x(w24_w,dpdn24_w);
    DF25_w= DFlux_dq_x(w25_w,dpdn25_w);
    DF26_w= DFlux_dq_x(w26_w,dpdn26_w);
    DF27_w= DFlux_dq_x(w27_w,dpdn27_w);
    DF01_n= DFlux_dq_y(w01_n,dpdn01_n);
    DF02_n= DFlux_dq_y(w02_n,dpdn02_n);
    DF03_n= DFlux_dq_y(w03_n,dpdn03_n);
    DF04_n= DFlux_dq_y(w04_n,dpdn04_n);
    DF05_n= DFlux_dq_y(w05_n,dpdn05_n);
    DF06_n= DFlux_dq_y(w06_n,dpdn06_n);
    DF07_n= DFlux_dq_y(w07_n,dpdn07_n);
    DF08_n= DFlux_dq_y(w08_n,dpdn08_n);
    DF09_n= DFlux_dq_y(w09_n,dpdn09_n);
    DF10_n= DFlux_dq_y(w10_n,dpdn10_n);
    DF11_n= DFlux_dq_y(w11_n,dpdn11_n);
    DF12_n= DFlux_dq_y(w12_n,dpdn12_n);
    DF13_n= DFlux_dq_y(w13_n,dpdn13_n);
    DF14_n= DFlux_dq_y(w14_n,dpdn14_n);
    DF15_n= DFlux_dq_y(w15_n,dpdn15_n);
    DF16_n= DFlux_dq_y(w16_n,dpdn16_n);
    DF17_n= DFlux_dq_y(w17_n,dpdn17_n);
    DF18_n= DFlux_dq_y(w18_n,dpdn18_n);
    DF19_n= DFlux_dq_y(w19_n,dpdn19_n);
    DF20_n= DFlux_dq_y(w20_n,dpdn20_n);
    DF21_n= DFlux_dq_y(w21_n,dpdn21_n);
    DF22_n= DFlux_dq_y(w22_n,dpdn22_n);
    DF23_n= DFlux_dq_y(w23_n,dpdn23_n);
    DF24_n= DFlux_dq_y(w24_n,dpdn24_n);
    DF25_n= DFlux_dq_y(w25_n,dpdn25_n);
    DF26_n= DFlux_dq_y(w26_n,dpdn26_n);
    DF27_n= DFlux_dq_y(w27_n,dpdn27_n);
    DF01_s= DFlux_dq_y(w01_s,dpdn01_s);
    DF02_s= DFlux_dq_y(w02_s,dpdn02_s);
    DF03_s= DFlux_dq_y(w03_s,dpdn03_s);
    DF04_s= DFlux_dq_y(w04_s,dpdn04_s);
    DF05_s= DFlux_dq_y(w05_s,dpdn05_s);
    DF06_s= DFlux_dq_y(w06_s,dpdn06_s);
    DF07_s= DFlux_dq_y(w07_s,dpdn07_s);
    DF08_s= DFlux_dq_y(w08_s,dpdn08_s);
    DF09_s= DFlux_dq_y(w09_s,dpdn09_s);
    DF10_s= DFlux_dq_y(w10_s,dpdn10_s);
    DF11_s= DFlux_dq_y(w11_s,dpdn11_s);
    DF12_s= DFlux_dq_y(w12_s,dpdn12_s);
    DF13_s= DFlux_dq_y(w13_s,dpdn13_s);
    DF14_s= DFlux_dq_y(w14_s,dpdn14_s);
    DF15_s= DFlux_dq_y(w15_s,dpdn15_s);
    DF16_s= DFlux_dq_y(w16_s,dpdn16_s);
    DF17_s= DFlux_dq_y(w17_s,dpdn17_s);
    DF18_s= DFlux_dq_y(w18_s,dpdn18_s);
    DF19_s= DFlux_dq_y(w19_s,dpdn19_s);
    DF20_s= DFlux_dq_y(w20_s,dpdn20_s);
    DF21_s= DFlux_dq_y(w21_s,dpdn21_s);
    DF22_s= DFlux_dq_y(w22_s,dpdn22_s);
    DF23_s= DFlux_dq_y(w23_s,dpdn23_s);
    DF24_s= DFlux_dq_y(w24_s,dpdn24_s);
    DF25_s= DFlux_dq_y(w25_s,dpdn25_s);
    DF26_s= DFlux_dq_y(w26_s,dpdn26_s);
    DF27_s= DFlux_dq_y(w27_s,dpdn27_s);
    DF01_u= DFlux_dq_z(w01_u,dpdn01_u);
    DF02_u= DFlux_dq_z(w02_u,dpdn02_u);
    DF03_u= DFlux_dq_z(w03_u,dpdn03_u);
    DF04_u= DFlux_dq_z(w04_u,dpdn04_u);
    DF05_u= DFlux_dq_z(w05_u,dpdn05_u);
    DF06_u= DFlux_dq_z(w06_u,dpdn06_u);
    DF07_u= DFlux_dq_z(w07_u,dpdn07_u);
    DF08_u= DFlux_dq_z(w08_u,dpdn08_u);
    DF09_u= DFlux_dq_z(w09_u,dpdn09_u);
    DF10_u= DFlux_dq_z(w10_u,dpdn10_u);
    DF11_u= DFlux_dq_z(w11_u,dpdn11_u);
    DF12_u= DFlux_dq_z(w12_u,dpdn12_u);
    DF13_u= DFlux_dq_z(w13_u,dpdn13_u);
    DF14_u= DFlux_dq_z(w14_u,dpdn14_u);
    DF15_u= DFlux_dq_z(w15_u,dpdn15_u);
    DF16_u= DFlux_dq_z(w16_u,dpdn16_u);
    DF17_u= DFlux_dq_z(w17_u,dpdn17_u);
    DF18_u= DFlux_dq_z(w18_u,dpdn18_u);
    DF19_u= DFlux_dq_z(w19_u,dpdn19_u);
    DF20_u= DFlux_dq_z(w20_u,dpdn20_u);
    DF21_u= DFlux_dq_z(w21_u,dpdn21_u);
    DF22_u= DFlux_dq_z(w22_u,dpdn22_u);
    DF23_u= DFlux_dq_z(w23_u,dpdn23_u);
    DF24_u= DFlux_dq_z(w24_u,dpdn24_u);
    DF25_u= DFlux_dq_z(w25_u,dpdn25_u);
    DF26_u= DFlux_dq_z(w26_u,dpdn26_u);
    DF27_u= DFlux_dq_z(w27_u,dpdn27_u);
    DF01_d= DFlux_dq_z(w01_d,dpdn01_d);
    DF02_d= DFlux_dq_z(w02_d,dpdn02_d);
    DF03_d= DFlux_dq_z(w03_d,dpdn03_d);
    DF04_d= DFlux_dq_z(w04_d,dpdn04_d);
    DF05_d= DFlux_dq_z(w05_d,dpdn05_d);
    DF06_d= DFlux_dq_z(w06_d,dpdn06_d);
    DF07_d= DFlux_dq_z(w07_d,dpdn07_d);
    DF08_d= DFlux_dq_z(w08_d,dpdn08_d);
    DF09_d= DFlux_dq_z(w09_d,dpdn09_d);
    DF10_d= DFlux_dq_z(w10_d,dpdn10_d);
    DF11_d= DFlux_dq_z(w11_d,dpdn11_d);
    DF12_d= DFlux_dq_z(w12_d,dpdn12_d);
    DF13_d= DFlux_dq_z(w13_d,dpdn13_d);
    DF14_d= DFlux_dq_z(w14_d,dpdn14_d);
    DF15_d= DFlux_dq_z(w15_d,dpdn15_d);
    DF16_d= DFlux_dq_z(w16_d,dpdn16_d);
    DF17_d= DFlux_dq_z(w17_d,dpdn17_d);
    DF18_d= DFlux_dq_z(w18_d,dpdn18_d);
    DF19_d= DFlux_dq_z(w19_d,dpdn19_d);
    DF20_d= DFlux_dq_z(w20_d,dpdn20_d);
    DF21_d= DFlux_dq_z(w21_d,dpdn21_d);
    DF22_d= DFlux_dq_z(w22_d,dpdn22_d);
    DF23_d= DFlux_dq_z(w23_d,dpdn23_d);
    DF24_d= DFlux_dq_z(w24_d,dpdn24_d);
    DF25_d= DFlux_dq_z(w25_d,dpdn25_d);
    DF26_d= DFlux_dq_z(w26_d,dpdn26_d);
    DF27_d= DFlux_dq_z(w27_d,dpdn27_d);

    east_trunc_311 = east_trunc_311 + psieT*(F03_e*wgtnuv1);
    east_trunc_321 = east_trunc_321 + psieT*(F06_e*wgtnuv1);
    east_trunc_331 = east_trunc_331 + psieT*(F09_e*wgtnuv1);
    east_trunc_312 = east_trunc_312 + psieT*(F12_e*wgtnuv1);
    east_trunc_322 = east_trunc_322 + psieT*(F15_e*wgtnuv1);
    east_trunc_332 = east_trunc_332 + psieT*(F18_e*wgtnuv1);
    east_trunc_313 = east_trunc_313 + psieT*(F21_e*wgtnuv1);
    east_trunc_323 = east_trunc_323 + psieT*(F24_e*wgtnuv1);
    east_trunc_333 = east_trunc_333 + psieT*(F27_e*wgtnuv1);        
    
    Jac_east_trunc_311 = Jac_east_trunc_311 + psieT*DF03_e*psie*wgtnuv1;
    Jac_east_trunc_321 = Jac_east_trunc_321 + psieT*DF06_e*psie*wgtnuv1;
    Jac_east_trunc_331 = Jac_east_trunc_331 + psieT*DF09_e*psie*wgtnuv1;          
    Jac_east_trunc_312 = Jac_east_trunc_312 + psieT*DF12_e*psie*wgtnuv1;
    Jac_east_trunc_322 = Jac_east_trunc_322 + psieT*DF15_e*psie*wgtnuv1;
    Jac_east_trunc_332 = Jac_east_trunc_332 + psieT*DF18_e*psie*wgtnuv1;      
    Jac_east_trunc_313 = Jac_east_trunc_313 + psieT*DF21_e*psie*wgtnuv1;
    Jac_east_trunc_323 = Jac_east_trunc_323 + psieT*DF24_e*psie*wgtnuv1;
    Jac_east_trunc_333 = Jac_east_trunc_333 + psieT*DF27_e*psie*wgtnuv1;        
   
       
    west_trunc_111 = west_trunc_111 - psiwT*(F01_w*wgtnuv1);
    west_trunc_121 = west_trunc_121 - psiwT*(F04_w*wgtnuv1);
    west_trunc_131 = west_trunc_131 - psiwT*(F07_w*wgtnuv1);           
    west_trunc_112 = west_trunc_112 - psiwT*(F10_w*wgtnuv1);
    west_trunc_122 = west_trunc_122 - psiwT*(F13_w*wgtnuv1);
    west_trunc_132 = west_trunc_132 - psiwT*(F16_w*wgtnuv1);           
    west_trunc_113 = west_trunc_113 - psiwT*(F19_w*wgtnuv1);
    west_trunc_123 = west_trunc_123 - psiwT*(F22_w*wgtnuv1);
    west_trunc_133 = west_trunc_133 - psiwT*(F25_w*wgtnuv1);
        
    Jac_west_trunc_111 = Jac_west_trunc_111 - psiwT*(DF01_w*psiw*wgtnuv1);
    Jac_west_trunc_121 = Jac_west_trunc_121 - psiwT*(DF04_w*psiw*wgtnuv1);
    Jac_west_trunc_131 = Jac_west_trunc_131 - psiwT*(DF07_w*psiw*wgtnuv1);         
    Jac_west_trunc_112 = Jac_west_trunc_112 - psiwT*(DF10_w*psiw*wgtnuv1);
    Jac_west_trunc_122 = Jac_west_trunc_122 - psiwT*(DF13_w*psiw*wgtnuv1);
    Jac_west_trunc_132 = Jac_west_trunc_132 - psiwT*(DF16_w*psiw*wgtnuv1);
    Jac_west_trunc_113 = Jac_west_trunc_113 - psiwT*(DF19_w*psiw*wgtnuv1);
    Jac_west_trunc_123 = Jac_west_trunc_123 - psiwT*(DF22_w*psiw*wgtnuv1);
    Jac_west_trunc_133 = Jac_west_trunc_133 - psiwT*(DF25_w*psiw*wgtnuv1);
    
    nort_trunc_131 = nort_trunc_131 + psinT*(F07_n*wgtnuv2);
    nort_trunc_231 = nort_trunc_231 + psinT*(F08_n*wgtnuv2);
    nort_trunc_331 = nort_trunc_331 + psinT*(F09_n*wgtnuv2);
    nort_trunc_132 = nort_trunc_132 + psinT*(F16_n*wgtnuv2);
    nort_trunc_232 = nort_trunc_232 + psinT*(F17_n*wgtnuv2);
    nort_trunc_332 = nort_trunc_332 + psinT*(F18_n*wgtnuv2);
    nort_trunc_133 = nort_trunc_133 + psinT*(F25_n*wgtnuv2);
    nort_trunc_233 = nort_trunc_233 + psinT*(F26_n*wgtnuv2);
    nort_trunc_333 = nort_trunc_333 + psinT*(F27_n*wgtnuv2);

    Jac_nort_trunc_131 = Jac_nort_trunc_131 + psinT*(DF07_n*psin*wgtnuv2);
    Jac_nort_trunc_231 = Jac_nort_trunc_231 + psinT*(DF08_n*psin*wgtnuv2);
    Jac_nort_trunc_331 = Jac_nort_trunc_331 + psinT*(DF09_n*psin*wgtnuv2);
    Jac_nort_trunc_132 = Jac_nort_trunc_132 + psinT*(DF16_n*psin*wgtnuv2);
    Jac_nort_trunc_232 = Jac_nort_trunc_232 + psinT*(DF17_n*psin*wgtnuv2);
    Jac_nort_trunc_332 = Jac_nort_trunc_332 + psinT*(DF18_n*psin*wgtnuv2);
    Jac_nort_trunc_133 = Jac_nort_trunc_133 + psinT*(DF25_n*psin*wgtnuv2);
    Jac_nort_trunc_233 = Jac_nort_trunc_233 + psinT*(DF26_n*psin*wgtnuv2);
    Jac_nort_trunc_333 = Jac_nort_trunc_333 + psinT*(DF27_n*psin*wgtnuv2);    
    
    sout_trunc_111 = sout_trunc_111 + psisT*(F01_s*wgtnuv2);
    sout_trunc_211 = sout_trunc_211 + psisT*(F02_s*wgtnuv2);
    sout_trunc_311 = sout_trunc_311 + psisT*(F03_s*wgtnuv2);
    sout_trunc_112 = sout_trunc_112 + psisT*(F10_s*wgtnuv2);
    sout_trunc_212 = sout_trunc_212 + psisT*(F11_s*wgtnuv2);
    sout_trunc_312 = sout_trunc_312 + psisT*(F12_s*wgtnuv2);
    sout_trunc_113 = sout_trunc_113 + psisT*(F19_s*wgtnuv2);
    sout_trunc_213 = sout_trunc_213 + psisT*(F20_s*wgtnuv2);
    sout_trunc_313 = sout_trunc_313 + psisT*(F21_s*wgtnuv2);

    Jac_sout_trunc_111 = Jac_sout_trunc_111 + psisT*(DF01_s*psis*wgtnuv2);
    Jac_sout_trunc_211 = Jac_sout_trunc_211 + psisT*(DF02_s*psis*wgtnuv2);
    Jac_sout_trunc_311 = Jac_sout_trunc_311 + psisT*(DF03_s*psis*wgtnuv2);
    Jac_sout_trunc_112 = Jac_sout_trunc_112 + psisT*(DF10_s*psis*wgtnuv2);
    Jac_sout_trunc_212 = Jac_sout_trunc_212 + psisT*(DF11_s*psis*wgtnuv2);
    Jac_sout_trunc_312 = Jac_sout_trunc_312 + psisT*(DF12_s*psis*wgtnuv2);
    Jac_sout_trunc_113 = Jac_sout_trunc_113 + psisT*(DF19_s*psis*wgtnuv2);
    Jac_sout_trunc_213 = Jac_sout_trunc_213 + psisT*(DF20_s*psis*wgtnuv2);
    Jac_sout_trunc_313 = Jac_sout_trunc_313 + psisT*(DF21_s*psis*wgtnuv2);    
        
    uppr_trunc_213  = uppr_trunc_213  + psiuT*(F20_u*wgtnuv3);
    uppr_trunc_313  = uppr_trunc_313  + psiuT*(F21_u*wgtnuv3);
    uppr_trunc_123  = uppr_trunc_123  + psiuT*(F22_u*wgtnuv3);
    uppr_trunc_223  = uppr_trunc_223  + psiuT*(F23_u*wgtnuv3);
    uppr_trunc_323  = uppr_trunc_323  + psiuT*(F24_u*wgtnuv3);
    uppr_trunc_133  = uppr_trunc_133  + psiuT*(F25_u*wgtnuv3);
    uppr_trunc_233  = uppr_trunc_233  + psiuT*(F26_u*wgtnuv3);
    uppr_trunc_333  = uppr_trunc_333  + psiuT*(F27_u*wgtnuv3);

    Jac_uppr_trunc_113 = Jac_uppr_trunc_113 + psiuT*(DF19_u*psiu*wgtnuv3);
    Jac_uppr_trunc_213 = Jac_uppr_trunc_213 + psiuT*(DF20_u*psiu*wgtnuv3);
    Jac_uppr_trunc_313 = Jac_uppr_trunc_313 + psiuT*(DF21_u*psiu*wgtnuv3);
    Jac_uppr_trunc_123 = Jac_uppr_trunc_123 + psiuT*(DF22_u*psiu*wgtnuv3);
    Jac_uppr_trunc_223 = Jac_uppr_trunc_223 + psiuT*(DF23_u*psiu*wgtnuv3);
    Jac_uppr_trunc_323 = Jac_uppr_trunc_323 + psiuT*(DF24_u*psiu*wgtnuv3);
    Jac_uppr_trunc_133 = Jac_uppr_trunc_133 + psiuT*(DF25_u*psiu*wgtnuv3);
    Jac_uppr_trunc_233 = Jac_uppr_trunc_233 + psiuT*(DF26_u*psiu*wgtnuv3);
    Jac_uppr_trunc_333 = Jac_uppr_trunc_333 + psiuT*(DF27_u*psiu*wgtnuv3);


    down_trunc_111  = down_trunc_111  + psidT*(F01_d*wgtnuv3);
    down_trunc_211  = down_trunc_211  + psidT*(F02_d*wgtnuv3);
    down_trunc_311  = down_trunc_311  + psidT*(F03_d*wgtnuv3);
    down_trunc_121  = down_trunc_121  + psidT*(F04_d*wgtnuv3);
    down_trunc_221  = down_trunc_221  + psidT*(F05_d*wgtnuv3);
    down_trunc_321  = down_trunc_321  + psidT*(F06_d*wgtnuv3);
    down_trunc_131  = down_trunc_131  + psidT*(F07_d*wgtnuv3);
    down_trunc_231  = down_trunc_231  + psidT*(F08_d*wgtnuv3);
    down_trunc_331  = down_trunc_331  + psidT*(F09_d*wgtnuv3);

    Jac_down_trunc_111 = Jac_down_trunc_111 + psidT*(DF01_d*psid*wgtnuv3);
    Jac_down_trunc_211 = Jac_down_trunc_211 + psidT*(DF02_d*psid*wgtnuv3);
    Jac_down_trunc_311 = Jac_down_trunc_311 + psidT*(DF03_d*psid*wgtnuv3);
    Jac_down_trunc_121 = Jac_down_trunc_121 + psidT*(DF04_d*psid*wgtnuv3);
    Jac_down_trunc_221 = Jac_down_trunc_221 + psidT*(DF05_d*psid*wgtnuv3);
    Jac_down_trunc_321 = Jac_down_trunc_321 + psidT*(DF06_d*psid*wgtnuv3);
    Jac_down_trunc_131 = Jac_down_trunc_131 + psidT*(DF07_d*psid*wgtnuv3);
    Jac_down_trunc_231 = Jac_down_trunc_231 + psidT*(DF08_d*psid*wgtnuv3);
    Jac_down_trunc_331 = Jac_down_trunc_331 + psidT*(DF09_d*psid*wgtnuv3);
    
    %Interface 1-2
    % cell left111
    %cell rite211
    %e1w
    %eastwest
    speedleft=speed_x(w01_e, dpdn01_e);
    speedavg=speed_x((w01_e+w02_w)/2, (dpdn01_e+dpdn02_w)/2 ) ;
    speedrite=speed_x(w02_w, dpdn02_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F01_e+F02_w-speed*(w02_w-w01_e))*wgtnuv1/2;
    F_dql=(DF01_e+speed)*wgtnuv1/2;
    F_dqr=(DF02_w-speed)*wgtnuv1/2;
    residual_cell_111=residual_cell_111+psieT*Flux;
    residual_cell_211=residual_cell_211-psiwT*Flux;
    Jac_cell_111=Jac_cell_111+psieT*F_dql*psie;
    Jac_east_othr_111=Jac_east_othr_111+psieT*F_dqr*psiw;
    Jac_cell_211=Jac_cell_211-psiwT*F_dqr*psiw;
    Jac_west_othr_211=Jac_west_othr_211-psiwT*F_dql*psie;

    %Interface 2-3
    % cell left211
    %cell rite311
    %e1w
    %eastwest
    speedleft=speed_x(w02_e, dpdn02_e);
    speedavg=speed_x((w02_e+w03_w)/2, (dpdn02_e+dpdn03_w)/2 ) ;
    speedrite=speed_x(w03_w, dpdn03_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F02_e+F03_w-speed*(w03_w-w02_e))*wgtnuv1/2;
    F_dql=(DF02_e+speed)*wgtnuv1/2;
    F_dqr=(DF03_w-speed)*wgtnuv1/2;
    residual_cell_211=residual_cell_211+psieT*Flux;
    residual_cell_311=residual_cell_311-psiwT*Flux;
    Jac_cell_211=Jac_cell_211+psieT*F_dql*psie;
    Jac_east_othr_211=Jac_east_othr_211+psieT*F_dqr*psiw;
    Jac_cell_311=Jac_cell_311-psiwT*F_dqr*psiw;
    Jac_west_othr_311=Jac_west_othr_311-psiwT*F_dql*psie;

    %Interface 4-5
    % cell left121
    %cell rite221
    %e1w
    %eastwest
    speedleft=speed_x(w04_e, dpdn04_e);
    speedavg=speed_x((w04_e+w05_w)/2, (dpdn04_e+dpdn05_w)/2 ) ;
    speedrite=speed_x(w05_w, dpdn05_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F04_e+F05_w-speed*(w05_w-w04_e))*wgtnuv1/2;
    F_dql=(DF04_e+speed)*wgtnuv1/2;
    F_dqr=(DF05_w-speed)*wgtnuv1/2;
    residual_cell_121=residual_cell_121+psieT*Flux;
    residual_cell_221=residual_cell_221-psiwT*Flux;
    Jac_cell_121=Jac_cell_121+psieT*F_dql*psie;
    Jac_east_othr_121=Jac_east_othr_121+psieT*F_dqr*psiw;
    Jac_cell_221=Jac_cell_221-psiwT*F_dqr*psiw;
    Jac_west_othr_221=Jac_west_othr_221-psiwT*F_dql*psie;

    %Interface 5-6
    % cell left221
    %cell rite321
    %e1w
    %eastwest
    speedleft=speed_x(w05_e, dpdn05_e);
    speedavg=speed_x((w05_e+w06_w)/2, (dpdn05_e+dpdn06_w)/2 ) ;
    speedrite=speed_x(w06_w, dpdn06_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F05_e+F06_w-speed*(w06_w-w05_e))*wgtnuv1/2;
    F_dql=(DF05_e+speed)*wgtnuv1/2;
    F_dqr=(DF06_w-speed)*wgtnuv1/2;
    residual_cell_221=residual_cell_221+psieT*Flux;
    residual_cell_321=residual_cell_321-psiwT*Flux;
    Jac_cell_221=Jac_cell_221+psieT*F_dql*psie;
    Jac_east_othr_221=Jac_east_othr_221+psieT*F_dqr*psiw;
    Jac_cell_321=Jac_cell_321-psiwT*F_dqr*psiw;
    Jac_west_othr_321=Jac_west_othr_321-psiwT*F_dql*psie;

    %Interface 7-8
    % cell left131
    %cell rite231
    %e1w
    %eastwest
    speedleft=speed_x(w07_e, dpdn07_e);
    speedavg=speed_x((w07_e+w08_w)/2, (dpdn07_e+dpdn08_w)/2 ) ;
    speedrite=speed_x(w08_w, dpdn08_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F07_e+F08_w-speed*(w08_w-w07_e))*wgtnuv1/2;
    F_dql=(DF07_e+speed)*wgtnuv1/2;
    F_dqr=(DF08_w-speed)*wgtnuv1/2;
    residual_cell_131=residual_cell_131+psieT*Flux;
    residual_cell_231=residual_cell_231-psiwT*Flux;
    Jac_cell_131=Jac_cell_131+psieT*F_dql*psie;
    Jac_east_othr_131=Jac_east_othr_131+psieT*F_dqr*psiw;
    Jac_cell_231=Jac_cell_231-psiwT*F_dqr*psiw;
    Jac_west_othr_231=Jac_west_othr_231-psiwT*F_dql*psie;

    %Interface 8-9
    % cell left231
    %cell rite331
    %e1w
    %eastwest
    speedleft=speed_x(w08_e, dpdn08_e);
    speedavg=speed_x((w08_e+w09_w)/2, (dpdn08_e+dpdn09_w)/2 ) ;
    speedrite=speed_x(w09_w, dpdn09_w);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F08_e+F09_w-speed*(w09_w-w08_e))*wgtnuv1/2;
    F_dql=(DF08_e+speed)*wgtnuv1/2;
    F_dqr=(DF09_w-speed)*wgtnuv1/2;
    residual_cell_231=residual_cell_231+psieT*Flux;
    residual_cell_331=residual_cell_331-psiwT*Flux;
    Jac_cell_231=Jac_cell_231+psieT*F_dql*psie;
    Jac_east_othr_231=Jac_east_othr_231+psieT*F_dqr*psiw;
    Jac_cell_331=Jac_cell_331-psiwT*F_dqr*psiw;
    Jac_west_othr_331=Jac_west_othr_331-psiwT*F_dql*psie;

    %Interface 10-11
    % cell left112
    %cell rite212
    %e1w
    %eastwest
    speedleft=speed_x(w10_e, dpdn10_e);
    speedavg=speed_x((w10_e+w11_w)/2, (dpdn10_e+dpdn11_w)/2 ) ;
    speedrite=speed_x(w11_w, dpdn11_w);
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

    %Interface 11-12
    % cell left212
    %cell rite312
    %e1w
    %eastwest
    speedleft=speed_x(w11_e, dpdn11_e);
    speedavg=speed_x((w11_e+w12_w)/2, (dpdn11_e+dpdn12_w)/2 ) ;
    speedrite=speed_x(w12_w, dpdn12_w);
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

    %Interface 13-14
    % cell left122
    %cell rite222
    %e1w
    %eastwest
    speedleft=speed_x(w13_e, dpdn13_e);
    speedavg=speed_x((w13_e+w14_w)/2, (dpdn13_e+dpdn14_w)/2 ) ;
    speedrite=speed_x(w14_w, dpdn14_w);
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

    %Interface 14-15
    % cell left222
    %cell rite322
    %e1w
    %eastwest
    speedleft=speed_x(w14_e, dpdn14_e);
    speedavg=speed_x((w14_e+w15_w)/2, (dpdn14_e+dpdn15_w)/2 ) ;
    speedrite=speed_x(w15_w, dpdn15_w);
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

    %Interface 16-17
    % cell left132
    %cell rite232
    %e1w
    %eastwest
    speedleft=speed_x(w16_e, dpdn16_e);
    speedavg=speed_x((w16_e+w17_w)/2, (dpdn16_e+dpdn17_w)/2 ) ;
    speedrite=speed_x(w17_w, dpdn17_w);
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

    %Interface 17-18
    % cell left232
    %cell rite332
    %e1w
    %eastwest
    speedleft=speed_x(w17_e, dpdn17_e);
    speedavg=speed_x((w17_e+w18_w)/2, (dpdn17_e+dpdn18_w)/2 ) ;
    speedrite=speed_x(w18_w, dpdn18_w);
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

    %Interface 19-20
    % cell left113
    %cell rite213
    %e1w
    %eastwest
    speedleft=speed_x(w19_e, dpdn19_e);
    speedavg=speed_x((w19_e+w20_w)/2, (dpdn19_e+dpdn20_w)/2 ) ;
    speedrite=speed_x(w20_w, dpdn20_w);
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

    %Interface 20-21
    % cell left213
    %cell rite313
    %e1w
    %eastwest
    speedleft=speed_x(w20_e, dpdn20_e);
    speedavg=speed_x((w20_e+w21_w)/2, (dpdn20_e+dpdn21_w)/2 ) ;
    speedrite=speed_x(w21_w, dpdn21_w);
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

    %Interface 22-23
    % cell left123
    %cell rite223
    %e1w
    %eastwest
    speedleft=speed_x(w22_e, dpdn22_e);
    speedavg=speed_x((w22_e+w23_w)/2, (dpdn22_e+dpdn23_w)/2 ) ;
    speedrite=speed_x(w23_w, dpdn23_w);
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

    %Interface 23-24
    % cell left223
    %cell rite323
    %e1w
    %eastwest
    speedleft=speed_x(w23_e, dpdn23_e);
    speedavg=speed_x((w23_e+w24_w)/2, (dpdn23_e+dpdn24_w)/2 ) ;
    speedrite=speed_x(w24_w, dpdn24_w);
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

    %Interface 25-26
    % cell left133
    %cell rite233
    %e1w
    %eastwest
    speedleft=speed_x(w25_e, dpdn25_e);
    speedavg=speed_x((w25_e+w26_w)/2, (dpdn25_e+dpdn26_w)/2 ) ;
    speedrite=speed_x(w26_w, dpdn26_w);
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

    %Interface 26-27
    % cell left233
    %cell rite333
    %e1w
    %eastwest
    speedleft=speed_x(w26_e, dpdn26_e);
    speedavg=speed_x((w26_e+w27_w)/2, (dpdn26_e+dpdn27_w)/2 ) ;
    speedrite=speed_x(w27_w, dpdn27_w);
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


    %Interface 1-4
    % cell left111
    %cell rite121
    %n2s
    %nortsout
    speedleft=speed_y(w01_n, dpdn01_n);
    speedavg=speed_y((w01_n+w04_s)/2, (dpdn01_n+dpdn04_s)/2 ) ;
    speedrite=speed_y(w04_s, dpdn04_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F01_n+F04_s-speed*(w04_s-w01_n))*wgtnuv2/2;
    F_dql=(DF01_n+speed)*wgtnuv2/2;
    F_dqr=(DF04_s-speed)*wgtnuv2/2;
    residual_cell_111=residual_cell_111+psinT*Flux;
    residual_cell_121=residual_cell_121-psisT*Flux;
    Jac_cell_111=Jac_cell_111+psinT*F_dql*psin;
    Jac_nort_othr_111=Jac_nort_othr_111+psinT*F_dqr*psis;
    Jac_cell_121=Jac_cell_121-psisT*F_dqr*psis;
    Jac_sout_othr_121=Jac_sout_othr_121-psisT*F_dql*psin;

    %Interface 4-7
    % cell left121
    %cell rite131
    %n2s
    %nortsout
    speedleft=speed_y(w04_n, dpdn04_n);
    speedavg=speed_y((w04_n+w07_s)/2, (dpdn04_n+dpdn07_s)/2 ) ;
    speedrite=speed_y(w07_s, dpdn07_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F04_n+F07_s-speed*(w07_s-w04_n))*wgtnuv2/2;
    F_dql=(DF04_n+speed)*wgtnuv2/2;
    F_dqr=(DF07_s-speed)*wgtnuv2/2;
    residual_cell_121=residual_cell_121+psinT*Flux;
    residual_cell_131=residual_cell_131-psisT*Flux;
    Jac_cell_121=Jac_cell_121+psinT*F_dql*psin;
    Jac_nort_othr_121=Jac_nort_othr_121+psinT*F_dqr*psis;
    Jac_cell_131=Jac_cell_131-psisT*F_dqr*psis;
    Jac_sout_othr_131=Jac_sout_othr_131-psisT*F_dql*psin;

    %Interface 2-5
    % cell left211
    %cell rite221
    %n2s
    %nortsout
    speedleft=speed_y(w02_n, dpdn02_n);
    speedavg=speed_y((w02_n+w05_s)/2, (dpdn02_n+dpdn05_s)/2 ) ;
    speedrite=speed_y(w05_s, dpdn05_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F02_n+F05_s-speed*(w05_s-w02_n))*wgtnuv2/2;
    F_dql=(DF02_n+speed)*wgtnuv2/2;
    F_dqr=(DF05_s-speed)*wgtnuv2/2;
    residual_cell_211=residual_cell_211+psinT*Flux;
    residual_cell_221=residual_cell_221-psisT*Flux;
    Jac_cell_211=Jac_cell_211+psinT*F_dql*psin;
    Jac_nort_othr_211=Jac_nort_othr_211+psinT*F_dqr*psis;
    Jac_cell_221=Jac_cell_221-psisT*F_dqr*psis;
    Jac_sout_othr_221=Jac_sout_othr_221-psisT*F_dql*psin;

    %Interface 5-8
    % cell left221
    %cell rite231
    %n2s
    %nortsout
    speedleft=speed_y(w05_n, dpdn05_n);
    speedavg=speed_y((w05_n+w08_s)/2, (dpdn05_n+dpdn08_s)/2 ) ;
    speedrite=speed_y(w08_s, dpdn08_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F05_n+F08_s-speed*(w08_s-w05_n))*wgtnuv2/2;
    F_dql=(DF05_n+speed)*wgtnuv2/2;
    F_dqr=(DF08_s-speed)*wgtnuv2/2;
    residual_cell_221=residual_cell_221+psinT*Flux;
    residual_cell_231=residual_cell_231-psisT*Flux;
    Jac_cell_221=Jac_cell_221+psinT*F_dql*psin;
    Jac_nort_othr_221=Jac_nort_othr_221+psinT*F_dqr*psis;
    Jac_cell_231=Jac_cell_231-psisT*F_dqr*psis;
    Jac_sout_othr_231=Jac_sout_othr_231-psisT*F_dql*psin;

    %Interface 3-6
    % cell left311
    %cell rite321
    %n2s
    %nortsout
    speedleft=speed_y(w03_n, dpdn03_n);
    speedavg=speed_y((w03_n+w06_s)/2, (dpdn03_n+dpdn06_s)/2 ) ;
    speedrite=speed_y(w06_s, dpdn06_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F03_n+F06_s-speed*(w06_s-w03_n))*wgtnuv2/2;
    F_dql=(DF03_n+speed)*wgtnuv2/2;
    F_dqr=(DF06_s-speed)*wgtnuv2/2;
    residual_cell_311=residual_cell_311+psinT*Flux;
    residual_cell_321=residual_cell_321-psisT*Flux;
    Jac_cell_311=Jac_cell_311+psinT*F_dql*psin;
    Jac_nort_othr_311=Jac_nort_othr_311+psinT*F_dqr*psis;
    Jac_cell_321=Jac_cell_321-psisT*F_dqr*psis;
    Jac_sout_othr_321=Jac_sout_othr_321-psisT*F_dql*psin;

    %Interface 6-9
    % cell left321
    %cell rite331
    %n2s
    %nortsout
    speedleft=speed_y(w06_n, dpdn06_n);
    speedavg=speed_y((w06_n+w09_s)/2, (dpdn06_n+dpdn09_s)/2 ) ;
    speedrite=speed_y(w09_s, dpdn09_s);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F06_n+F09_s-speed*(w09_s-w06_n))*wgtnuv2/2;
    F_dql=(DF06_n+speed)*wgtnuv2/2;
    F_dqr=(DF09_s-speed)*wgtnuv2/2;
    residual_cell_321=residual_cell_321+psinT*Flux;
    residual_cell_331=residual_cell_331-psisT*Flux;
    Jac_cell_321=Jac_cell_321+psinT*F_dql*psin;
    Jac_nort_othr_321=Jac_nort_othr_321+psinT*F_dqr*psis;
    Jac_cell_331=Jac_cell_331-psisT*F_dqr*psis;
    Jac_sout_othr_331=Jac_sout_othr_331-psisT*F_dql*psin;

    %Interface 10-13
    % cell left112
    %cell rite122
    %n2s
    %nortsout
    speedleft=speed_y(w10_n, dpdn10_n);
    speedavg=speed_y((w10_n+w13_s)/2, (dpdn10_n+dpdn13_s)/2 ) ;
    speedrite=speed_y(w13_s, dpdn13_s);
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

    %Interface 13-16
    % cell left122
    %cell rite132
    %n2s
    %nortsout
    speedleft=speed_y(w13_n, dpdn13_n);
    speedavg=speed_y((w13_n+w16_s)/2, (dpdn13_n+dpdn16_s)/2 ) ;
    speedrite=speed_y(w16_s, dpdn16_s);
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

    %Interface 11-14
    % cell left212
    %cell rite222
    %n2s
    %nortsout
    speedleft=speed_y(w11_n, dpdn11_n);
    speedavg=speed_y((w11_n+w14_s)/2, (dpdn11_n+dpdn14_s)/2 ) ;
    speedrite=speed_y(w14_s, dpdn14_s);
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

    %Interface 14-17
    % cell left222
    %cell rite232
    %n2s
    %nortsout
    speedleft=speed_y(w14_n, dpdn14_n);
    speedavg=speed_y((w14_n+w17_s)/2, (dpdn14_n+dpdn17_s)/2 ) ;
    speedrite=speed_y(w17_s, dpdn17_s);
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

    %Interface 12-15
    % cell left312
    %cell rite322
    %n2s
    %nortsout
    speedleft=speed_y(w12_n, dpdn12_n);
    speedavg=speed_y((w12_n+w15_s)/2, (dpdn12_n+dpdn15_s)/2 ) ;
    speedrite=speed_y(w15_s, dpdn15_s);
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

    %Interface 15-18
    % cell left322
    %cell rite332
    %n2s
    %nortsout
    speedleft=speed_y(w15_n, dpdn15_n);
    speedavg=speed_y((w15_n+w18_s)/2, (dpdn15_n+dpdn18_s)/2 ) ;
    speedrite=speed_y(w18_s, dpdn18_s);
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

    %Interface 19-22
    % cell left113
    %cell rite123
    %n2s
    %nortsout
    speedleft=speed_y(w19_n, dpdn19_n);
    speedavg=speed_y((w19_n+w22_s)/2, (dpdn19_n+dpdn22_s)/2 ) ;
    speedrite=speed_y(w22_s, dpdn22_s);
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

    %Interface 22-25
    % cell left123
    %cell rite133
    %n2s
    %nortsout
    speedleft=speed_y(w22_n, dpdn22_n);
    speedavg=speed_y((w22_n+w25_s)/2, (dpdn22_n+dpdn25_s)/2 ) ;
    speedrite=speed_y(w25_s, dpdn25_s);
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

    %Interface 20-23
    % cell left213
    %cell rite223
    %n2s
    %nortsout
    speedleft=speed_y(w20_n, dpdn20_n);
    speedavg=speed_y((w20_n+w23_s)/2, (dpdn20_n+dpdn23_s)/2 ) ;
    speedrite=speed_y(w23_s, dpdn23_s);
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

    %Interface 23-26
    % cell left223
    %cell rite233
    %n2s
    %nortsout
    speedleft=speed_y(w23_n, dpdn23_n);
    speedavg=speed_y((w23_n+w26_s)/2, (dpdn23_n+dpdn26_s)/2 ) ;
    speedrite=speed_y(w26_s, dpdn26_s);
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

    %Interface 21-24
    % cell left313
    %cell rite323
    %n2s
    %nortsout
    speedleft=speed_y(w21_n, dpdn21_n);
    speedavg=speed_y((w21_n+w24_s)/2, (dpdn21_n+dpdn24_s)/2 ) ;
    speedrite=speed_y(w24_s, dpdn24_s);
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

    %Interface 24-27
    % cell left323
    %cell rite333
    %n2s
    %nortsout
    speedleft=speed_y(w24_n, dpdn24_n);
    speedavg=speed_y((w24_n+w27_s)/2, (dpdn24_n+dpdn27_s)/2 ) ;
    speedrite=speed_y(w27_s, dpdn27_s);
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



    %Interface 1-10
    % cell left111
    %cell rite112
    %u3d
    %upprdown
    speedleft=speed_z(w01_u, dpdn01_u);
    speedavg=speed_z((w01_u+w10_d)/2, (dpdn01_u+dpdn10_d)/2 ) ;
    speedrite=speed_z(w10_d, dpdn10_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F01_u+F10_d-speed*(w10_d-w01_u))*wgtnuv3/2;
    F_dql=(DF01_u+speed)*wgtnuv3/2;
    F_dqr=(DF10_d-speed)*wgtnuv3/2;
    residual_cell_111=residual_cell_111+psiuT*Flux;
    residual_cell_112=residual_cell_112-psidT*Flux;
    Jac_cell_111=Jac_cell_111+psiuT*F_dql*psiu;
    Jac_uppr_othr_111=Jac_uppr_othr_111+psiuT*F_dqr*psid;
    Jac_cell_112=Jac_cell_112-psidT*F_dqr*psid;
    Jac_down_othr_112=Jac_down_othr_112-psidT*F_dql*psiu;

    %Interface 2-11
    % cell left211
    %cell rite212
    %u3d
    %upprdown
    speedleft=speed_z(w02_u, dpdn02_u);
    speedavg=speed_z((w02_u+w11_d)/2, (dpdn02_u+dpdn11_d)/2 ) ;
    speedrite=speed_z(w11_d, dpdn11_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F02_u+F11_d-speed*(w11_d-w02_u))*wgtnuv3/2;
    F_dql=(DF02_u+speed)*wgtnuv3/2;
    F_dqr=(DF11_d-speed)*wgtnuv3/2;
    residual_cell_211=residual_cell_211+psiuT*Flux;
    residual_cell_212=residual_cell_212-psidT*Flux;
    Jac_cell_211=Jac_cell_211+psiuT*F_dql*psiu;
    Jac_uppr_othr_211=Jac_uppr_othr_211+psiuT*F_dqr*psid;
    Jac_cell_212=Jac_cell_212-psidT*F_dqr*psid;
    Jac_down_othr_212=Jac_down_othr_212-psidT*F_dql*psiu;

    %Interface 3-12
    % cell left311
    %cell rite312
    %u3d
    %upprdown
    speedleft=speed_z(w03_u, dpdn03_u);
    speedavg=speed_z((w03_u+w12_d)/2, (dpdn03_u+dpdn12_d)/2 ) ;
    speedrite=speed_z(w12_d, dpdn12_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F03_u+F12_d-speed*(w12_d-w03_u))*wgtnuv3/2;
    F_dql=(DF03_u+speed)*wgtnuv3/2;
    F_dqr=(DF12_d-speed)*wgtnuv3/2;
    residual_cell_311=residual_cell_311+psiuT*Flux;
    residual_cell_312=residual_cell_312-psidT*Flux;
    Jac_cell_311=Jac_cell_311+psiuT*F_dql*psiu;
    Jac_uppr_othr_311=Jac_uppr_othr_311+psiuT*F_dqr*psid;
    Jac_cell_312=Jac_cell_312-psidT*F_dqr*psid;
    Jac_down_othr_312=Jac_down_othr_312-psidT*F_dql*psiu;

    %Interface 4-13
    % cell left121
    %cell rite122
    %u3d
    %upprdown
    speedleft=speed_z(w04_u, dpdn04_u);
    speedavg=speed_z((w04_u+w13_d)/2, (dpdn04_u+dpdn13_d)/2 ) ;
    speedrite=speed_z(w13_d, dpdn13_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F04_u+F13_d-speed*(w13_d-w04_u))*wgtnuv3/2;
    F_dql=(DF04_u+speed)*wgtnuv3/2;
    F_dqr=(DF13_d-speed)*wgtnuv3/2;
    residual_cell_121=residual_cell_121+psiuT*Flux;
    residual_cell_122=residual_cell_122-psidT*Flux;
    Jac_cell_121=Jac_cell_121+psiuT*F_dql*psiu;
    Jac_uppr_othr_121=Jac_uppr_othr_121+psiuT*F_dqr*psid;
    Jac_cell_122=Jac_cell_122-psidT*F_dqr*psid;
    Jac_down_othr_122=Jac_down_othr_122-psidT*F_dql*psiu;

    %Interface 5-14
    % cell left221
    %cell rite222
    %u3d
    %upprdown
    speedleft=speed_z(w05_u, dpdn05_u);
    speedavg=speed_z((w05_u+w14_d)/2, (dpdn05_u+dpdn14_d)/2 ) ;
    speedrite=speed_z(w14_d, dpdn14_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F05_u+F14_d-speed*(w14_d-w05_u))*wgtnuv3/2;
    F_dql=(DF05_u+speed)*wgtnuv3/2;
    F_dqr=(DF14_d-speed)*wgtnuv3/2;
    residual_cell_221=residual_cell_221+psiuT*Flux;
    residual_cell_222=residual_cell_222-psidT*Flux;
    Jac_cell_221=Jac_cell_221+psiuT*F_dql*psiu;
    Jac_uppr_othr_221=Jac_uppr_othr_221+psiuT*F_dqr*psid;
    Jac_cell_222=Jac_cell_222-psidT*F_dqr*psid;
    Jac_down_othr_222=Jac_down_othr_222-psidT*F_dql*psiu;

    %Interface 6-15
    % cell left321
    %cell rite322
    %u3d
    %upprdown
    speedleft=speed_z(w06_u, dpdn06_u);
    speedavg=speed_z((w06_u+w15_d)/2, (dpdn06_u+dpdn15_d)/2 ) ;
    speedrite=speed_z(w15_d, dpdn15_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F06_u+F15_d-speed*(w15_d-w06_u))*wgtnuv3/2;
    F_dql=(DF06_u+speed)*wgtnuv3/2;
    F_dqr=(DF15_d-speed)*wgtnuv3/2;
    residual_cell_321=residual_cell_321+psiuT*Flux;
    residual_cell_322=residual_cell_322-psidT*Flux;
    Jac_cell_321=Jac_cell_321+psiuT*F_dql*psiu;
    Jac_uppr_othr_321=Jac_uppr_othr_321+psiuT*F_dqr*psid;
    Jac_cell_322=Jac_cell_322-psidT*F_dqr*psid;
    Jac_down_othr_322=Jac_down_othr_322-psidT*F_dql*psiu;

    %Interface 7-16
    % cell left131
    %cell rite132
    %u3d
    %upprdown
    speedleft=speed_z(w07_u, dpdn07_u);
    speedavg=speed_z((w07_u+w16_d)/2, (dpdn07_u+dpdn16_d)/2 ) ;
    speedrite=speed_z(w16_d, dpdn16_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F07_u+F16_d-speed*(w16_d-w07_u))*wgtnuv3/2;
    F_dql=(DF07_u+speed)*wgtnuv3/2;
    F_dqr=(DF16_d-speed)*wgtnuv3/2;
    residual_cell_131=residual_cell_131+psiuT*Flux;
    residual_cell_132=residual_cell_132-psidT*Flux;
    Jac_cell_131=Jac_cell_131+psiuT*F_dql*psiu;
    Jac_uppr_othr_131=Jac_uppr_othr_131+psiuT*F_dqr*psid;
    Jac_cell_132=Jac_cell_132-psidT*F_dqr*psid;
    Jac_down_othr_132=Jac_down_othr_132-psidT*F_dql*psiu;

    %Interface 8-17
    % cell left231
    %cell rite232
    %u3d
    %upprdown
    speedleft=speed_z(w08_u, dpdn08_u);
    speedavg=speed_z((w08_u+w17_d)/2, (dpdn08_u+dpdn17_d)/2 ) ;
    speedrite=speed_z(w17_d, dpdn17_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F08_u+F17_d-speed*(w17_d-w08_u))*wgtnuv3/2;
    F_dql=(DF08_u+speed)*wgtnuv3/2;
    F_dqr=(DF17_d-speed)*wgtnuv3/2;
    residual_cell_231=residual_cell_231+psiuT*Flux;
    residual_cell_232=residual_cell_232-psidT*Flux;
    Jac_cell_231=Jac_cell_231+psiuT*F_dql*psiu;
    Jac_uppr_othr_231=Jac_uppr_othr_231+psiuT*F_dqr*psid;
    Jac_cell_232=Jac_cell_232-psidT*F_dqr*psid;
    Jac_down_othr_232=Jac_down_othr_232-psidT*F_dql*psiu;

    %Interface 9-18
    % cell left331
    %cell rite332
    %u3d
    %upprdown
    speedleft=speed_z(w09_u, dpdn09_u);
    speedavg=speed_z((w09_u+w18_d)/2, (dpdn09_u+dpdn18_d)/2 ) ;
    speedrite=speed_z(w18_d, dpdn18_d);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);
    Flux=(F09_u+F18_d-speed*(w18_d-w09_u))*wgtnuv3/2;
    F_dql=(DF09_u+speed)*wgtnuv3/2;
    F_dqr=(DF18_d-speed)*wgtnuv3/2;
    residual_cell_331=residual_cell_331+psiuT*Flux;
    residual_cell_332=residual_cell_332-psidT*Flux;
    Jac_cell_331=Jac_cell_331+psiuT*F_dql*psiu;
    Jac_uppr_othr_331=Jac_uppr_othr_331+psiuT*F_dqr*psid;
    Jac_cell_332=Jac_cell_332-psidT*F_dqr*psid;
    Jac_down_othr_332=Jac_down_othr_332-psidT*F_dql*psiu;

    %Interface 10-19
    % cell left112
    %cell rite113
    %u3d
    %upprdown
    speedleft=speed_z(w10_u, dpdn10_u);
    speedavg=speed_z((w10_u+w19_d)/2, (dpdn10_u+dpdn19_d)/2 ) ;
    speedrite=speed_z(w19_d, dpdn19_d);
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

    %Interface 11-20
    % cell left212
    %cell rite213
    %u3d
    %upprdown
    speedleft=speed_z(w11_u, dpdn11_u);
    speedavg=speed_z((w11_u+w20_d)/2, (dpdn11_u+dpdn20_d)/2 ) ;
    speedrite=speed_z(w20_d, dpdn20_d);
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

    %Interface 12-21
    % cell left312
    %cell rite313
    %u3d
    %upprdown
    speedleft=speed_z(w12_u, dpdn12_u);
    speedavg=speed_z((w12_u+w21_d)/2, (dpdn12_u+dpdn21_d)/2 ) ;
    speedrite=speed_z(w21_d, dpdn21_d);
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

    %Interface 13-22
    % cell left122
    %cell rite123
    %u3d
    %upprdown
    speedleft=speed_z(w13_u, dpdn13_u);
    speedavg=speed_z((w13_u+w22_d)/2, (dpdn13_u+dpdn22_d)/2 ) ;
    speedrite=speed_z(w22_d, dpdn22_d);
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

    %Interface 14-23
    % cell left222
    %cell rite223
    %u3d
    %upprdown
    speedleft=speed_z(w14_u, dpdn14_u);
    speedavg=speed_z((w14_u+w23_d)/2, (dpdn14_u+dpdn23_d)/2 ) ;
    speedrite=speed_z(w23_d, dpdn23_d);
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

    %Interface 15-24
    % cell left322
    %cell rite323
    %u3d
    %upprdown
    speedleft=speed_z(w15_u, dpdn15_u);
    speedavg=speed_z((w15_u+w24_d)/2, (dpdn15_u+dpdn24_d)/2 ) ;
    speedrite=speed_z(w24_d, dpdn24_d);
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

    %Interface 16-25
    % cell left132
    %cell rite133
    %u3d
    %upprdown
    speedleft=speed_z(w16_u, dpdn16_u);
    speedavg=speed_z((w16_u+w25_d)/2, (dpdn16_u+dpdn25_d)/2 ) ;
    speedrite=speed_z(w25_d, dpdn25_d);
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

    %Interface 17-26
    % cell left232
    %cell rite233
    %u3d
    %upprdown
    speedleft=speed_z(w17_u, dpdn17_u);
    speedavg=speed_z((w17_u+w26_d)/2, (dpdn17_u+dpdn26_d)/2 ) ;
    speedrite=speed_z(w26_d, dpdn26_d);
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

    %Interface 18-27
    % cell left332
    %cell rite333
    %u3d
    %upprdown
    speedleft=speed_z(w18_u, dpdn18_u);
    speedavg=speed_z((w18_u+w27_d)/2, (dpdn18_u+dpdn27_d)/2 ) ;
    speedrite=speed_z(w27_d, dpdn27_d);
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

Jacobian(indices(1),indices(1))=Jac_cell_111+Jac_east_flux_111+Jac_west_trunc_111+Jac_nort_flux_111+Jac_sout_trunc_111+Jac_uppr_flux_111+Jac_down_trunc_111;%
Jacobian(indices(2),indices(2))=Jac_cell_211+Jac_east_flux_211+Jac_west_flux_211+Jac_nort_flux_211+Jac_sout_trunc_211+Jac_uppr_flux_211+Jac_down_trunc_211;%
Jacobian(indices(3),indices(3))=Jac_cell_311+Jac_east_trunc_311+Jac_west_flux_311+Jac_nort_flux_311+Jac_sout_trunc_311+Jac_uppr_flux_311+Jac_down_trunc_311;%
Jacobian(indices(4),indices(4))=Jac_cell_121+Jac_east_flux_121+Jac_west_trunc_121+Jac_nort_flux_121+Jac_sout_flux_121+Jac_uppr_flux_121+Jac_down_trunc_121;%
Jacobian(indices(5),indices(5))=Jac_cell_221+Jac_east_flux_221+Jac_west_flux_221+Jac_nort_flux_221+Jac_sout_flux_221+Jac_uppr_flux_221+Jac_down_trunc_221;%
Jacobian(indices(6),indices(6))=Jac_cell_321+Jac_east_trunc_321+Jac_west_flux_321+Jac_nort_flux_321+Jac_sout_flux_321+Jac_uppr_flux_321+Jac_down_trunc_321;%
Jacobian(indices(7),indices(7))=Jac_cell_131+Jac_east_flux_131+Jac_west_trunc_131+Jac_nort_trunc_131+Jac_sout_flux_131+Jac_uppr_flux_131+Jac_down_trunc_131;%
Jacobian(indices(8),indices(8))=Jac_cell_231+Jac_east_flux_231+Jac_west_flux_231+Jac_nort_trunc_231+Jac_sout_flux_231+Jac_uppr_flux_231+Jac_down_trunc_231;%
Jacobian(indices(9),indices(9))=Jac_cell_331+Jac_east_trunc_331+Jac_west_flux_331+Jac_nort_trunc_331+Jac_sout_flux_331+Jac_uppr_flux_331+Jac_down_trunc_331;%
Jacobian(indices(10),indices(10))=Jac_cell_112+Jac_east_flux_112+Jac_west_trunc_112+Jac_nort_flux_112+Jac_sout_trunc_112+Jac_uppr_flux_112+Jac_down_flux_112;%
Jacobian(indices(11),indices(11))=Jac_cell_212+Jac_east_flux_212+Jac_west_flux_212+Jac_nort_flux_212+Jac_sout_trunc_212+Jac_uppr_flux_212+Jac_down_flux_212;%
Jacobian(indices(12),indices(12))=Jac_cell_312+Jac_east_trunc_312+Jac_west_flux_312+Jac_nort_flux_312+Jac_sout_trunc_312+Jac_uppr_flux_312+Jac_down_flux_312;%
Jacobian(indices(13),indices(13))=Jac_cell_122+Jac_east_flux_122+Jac_west_trunc_122+Jac_nort_flux_122+Jac_sout_flux_122+Jac_uppr_flux_122+Jac_down_flux_122;%
Jacobian(indices(14),indices(14))=Jac_cell_222+Jac_east_flux_222+Jac_west_flux_222+Jac_nort_flux_222+Jac_sout_flux_222+Jac_uppr_flux_222+Jac_down_flux_222;%
Jacobian(indices(15),indices(15))=Jac_cell_322+Jac_east_trunc_322+Jac_west_flux_322+Jac_nort_flux_322+Jac_sout_flux_322+Jac_uppr_flux_322+Jac_down_flux_322;%
Jacobian(indices(16),indices(16))=Jac_cell_132+Jac_east_flux_132+Jac_west_trunc_132+Jac_nort_trunc_132+Jac_sout_flux_132+Jac_uppr_flux_132+Jac_down_flux_132;%
Jacobian(indices(17),indices(17))=Jac_cell_232+Jac_east_flux_232+Jac_west_flux_232+Jac_nort_trunc_232+Jac_sout_flux_232+Jac_uppr_flux_232+Jac_down_flux_232;%
Jacobian(indices(18),indices(18))=Jac_cell_332+Jac_east_trunc_332+Jac_west_flux_332+Jac_nort_trunc_332+Jac_sout_flux_332+Jac_uppr_flux_332+Jac_down_flux_332;%
Jacobian(indices(19),indices(19))=Jac_cell_113+Jac_east_flux_113+Jac_west_trunc_113+Jac_nort_flux_113+Jac_sout_trunc_113+Jac_uppr_trunc_113+Jac_down_flux_113;%
Jacobian(indices(20),indices(20))=Jac_cell_213+Jac_east_flux_213+Jac_west_flux_213+Jac_nort_flux_213+Jac_sout_trunc_213+Jac_uppr_trunc_213+Jac_down_flux_213;%
Jacobian(indices(21),indices(21))=Jac_cell_313+Jac_east_trunc_313+Jac_west_flux_313+Jac_nort_flux_313+Jac_sout_trunc_313+Jac_uppr_trunc_313+Jac_down_flux_313;%
Jacobian(indices(22),indices(22))=Jac_cell_123+Jac_east_flux_123+Jac_west_trunc_123+Jac_nort_flux_123+Jac_sout_flux_123+Jac_uppr_trunc_123+Jac_down_flux_123;%
Jacobian(indices(23),indices(23))=Jac_cell_223+Jac_east_flux_223+Jac_west_flux_223+Jac_nort_flux_223+Jac_sout_flux_223+Jac_uppr_trunc_223+Jac_down_flux_223;%
Jacobian(indices(24),indices(24))=Jac_cell_323+Jac_east_trunc_323+Jac_west_flux_323+Jac_nort_flux_323+Jac_sout_flux_323+Jac_uppr_trunc_323+Jac_down_flux_323;%
Jacobian(indices(25),indices(25))=Jac_cell_133+Jac_east_flux_133+Jac_west_trunc_133+Jac_nort_trunc_133+Jac_sout_flux_133+Jac_uppr_trunc_133+Jac_down_flux_133;%
Jacobian(indices(26),indices(26))=Jac_cell_233+Jac_east_flux_233+Jac_west_flux_233+Jac_nort_trunc_233+Jac_sout_flux_233+Jac_uppr_trunc_233+Jac_down_flux_233;%
Jacobian(indices(27),indices(27))=Jac_cell_333+Jac_east_trunc_333+Jac_west_flux_333+Jac_nort_trunc_333+Jac_sout_flux_333+Jac_uppr_trunc_333+Jac_down_flux_333;%

Jacobian(indices(1),indices(2)) = Jac_east_othr_111;%
Jacobian(indices(2),indices(3)) = Jac_east_othr_211;%
%Jacobian(indices(3),indices(4)) = Jac_east_othr_311;%
Jacobian(indices(4),indices(5)) = Jac_east_othr_121;%
Jacobian(indices(5),indices(6)) = Jac_east_othr_221;%
%Jacobian(indices(6),indices(7)) = Jac_east_othr_321;%
Jacobian(indices(7),indices(8)) = Jac_east_othr_131;%
Jacobian(indices(8),indices(9)) = Jac_east_othr_231;%
%Jacobian(indices(9),indices(10)) = Jac_east_othr_331;%
Jacobian(indices(10),indices(11)) = Jac_east_othr_112;%
Jacobian(indices(11),indices(12)) = Jac_east_othr_212;%
%Jacobian(indices(12),indices(13)) = Jac_east_othr_312;%
Jacobian(indices(13),indices(14)) = Jac_east_othr_122;%
Jacobian(indices(14),indices(15)) = Jac_east_othr_222;%
%Jacobian(indices(15),indices(16)) = Jac_east_othr_322;%
Jacobian(indices(16),indices(17)) = Jac_east_othr_132;%
Jacobian(indices(17),indices(18)) = Jac_east_othr_232;%
%Jacobian(indices(18),indices(19)) = Jac_east_othr_332;%
Jacobian(indices(19),indices(20)) = Jac_east_othr_113;%
Jacobian(indices(20),indices(21)) = Jac_east_othr_213;%
%Jacobian(indices(21),indices(22)) = Jac_east_othr_313;%
Jacobian(indices(22),indices(23)) = Jac_east_othr_123;%
Jacobian(indices(23),indices(24)) = Jac_east_othr_223;%
%Jacobian(indices(24),indices(25)) = Jac_east_othr_323;%
Jacobian(indices(25),indices(26)) = Jac_east_othr_133;%
Jacobian(indices(26),indices(27)) = Jac_east_othr_233;%
%;%
%Jacobian(indices(1),indices(0)) = Jac_west_othr_111;%
Jacobian(indices(2),indices(1)) = Jac_west_othr_211;%
Jacobian(indices(3),indices(2)) = Jac_west_othr_311;%
%Jacobian(indices(4),indices(3)) = Jac_west_othr_121;%
Jacobian(indices(5),indices(4)) = Jac_west_othr_221;%
Jacobian(indices(6),indices(5)) = Jac_west_othr_321;%
%Jacobian(indices(7),indices(6)) = Jac_west_othr_131;%
Jacobian(indices(8),indices(7)) = Jac_west_othr_231;%
Jacobian(indices(9),indices(8)) = Jac_west_othr_331;%
%Jacobian(indices(10),indices(9)) = Jac_west_othr_112;%
Jacobian(indices(11),indices(10)) = Jac_west_othr_212;%
Jacobian(indices(12),indices(11)) = Jac_west_othr_312;%
%Jacobian(indices(13),indices(12)) = Jac_west_othr_122;%
Jacobian(indices(14),indices(13)) = Jac_west_othr_222;%
Jacobian(indices(15),indices(14)) = Jac_west_othr_322;%
%Jacobian(indices(16),indices(15)) = Jac_west_othr_132;%
Jacobian(indices(17),indices(16)) = Jac_west_othr_232;%
Jacobian(indices(18),indices(17)) = Jac_west_othr_332;%
%Jacobian(indices(19),indices(18)) = Jac_west_othr_113;%
Jacobian(indices(20),indices(19)) = Jac_west_othr_213;%
Jacobian(indices(21),indices(20)) = Jac_west_othr_313;%
%Jacobian(indices(22),indices(21)) = Jac_west_othr_123;%
Jacobian(indices(23),indices(22)) = Jac_west_othr_223;%
Jacobian(indices(24),indices(23)) = Jac_west_othr_323;%
%Jacobian(indices(25),indices(24)) = Jac_west_othr_133;%
Jacobian(indices(26),indices(25)) = Jac_west_othr_233;%
Jacobian(indices(27),indices(26)) = Jac_west_othr_333;%
%;%
Jacobian(indices(1),indices(4)) = Jac_nort_othr_111;%
Jacobian(indices(2),indices(5)) = Jac_nort_othr_211;%
Jacobian(indices(3),indices(6)) = Jac_nort_othr_311;%
Jacobian(indices(4),indices(7)) = Jac_nort_othr_121;%
Jacobian(indices(5),indices(8)) = Jac_nort_othr_221;%
Jacobian(indices(6),indices(9)) = Jac_nort_othr_321;%
%Jacobian(indices(7),indices(10)) = Jac_nort_othr_131;%
%Jacobian(indices(8),indices(11)) = Jac_nort_othr_231;%
%Jacobian(indices(9),indices(12)) = Jac_nort_othr_331;%
Jacobian(indices(10),indices(13)) = Jac_nort_othr_112;%
Jacobian(indices(11),indices(14)) = Jac_nort_othr_212;%
Jacobian(indices(12),indices(15)) = Jac_nort_othr_312;%
Jacobian(indices(13),indices(16)) = Jac_nort_othr_122;%
Jacobian(indices(14),indices(17)) = Jac_nort_othr_222;%
Jacobian(indices(15),indices(18)) = Jac_nort_othr_322;%
%Jacobian(indices(16),indices(19)) = Jac_nort_othr_132;%
%Jacobian(indices(17),indices(20)) = Jac_nort_othr_232;%
%Jacobian(indices(18),indices(21)) = Jac_nort_othr_332;%
Jacobian(indices(19),indices(22)) = Jac_nort_othr_113;%
Jacobian(indices(20),indices(23)) = Jac_nort_othr_213;%
Jacobian(indices(21),indices(24)) = Jac_nort_othr_313;%
Jacobian(indices(22),indices(25)) = Jac_nort_othr_123;%
Jacobian(indices(23),indices(26)) = Jac_nort_othr_223;%
Jacobian(indices(24),indices(27)) = Jac_nort_othr_323;%
%Jacobian(indices(25),indices(28)) = Jac_nort_othr_133;%
%Jacobian(indices(26),indices(29)) = Jac_nort_othr_233;%
%;%
%Jacobian(indices(1),indices(-15)) = Jac_sout_othr_111;%
%Jacobian(indices(2),indices(-1)) = Jac_sout_othr_211;%
%Jacobian(indices(3),indices(0)) = Jac_sout_othr_311;%
%Jacobian(indices(4),indices(1)) = Jac_sout_othr_121;%
Jacobian(indices(5),indices(2)) = Jac_sout_othr_221;%
Jacobian(indices(6),indices(3)) = Jac_sout_othr_321;%
Jacobian(indices(7),indices(4)) = Jac_sout_othr_131;%
Jacobian(indices(8),indices(5)) = Jac_sout_othr_231;%
Jacobian(indices(9),indices(6)) = Jac_sout_othr_331;%
Jacobian(indices(10),indices(7)) = Jac_sout_othr_112;%
%Jacobian(indices(11),indices(8)) = Jac_sout_othr_212;%
%Jacobian(indices(12),indices(9)) = Jac_sout_othr_312;%
%Jacobian(indices(13),indices(10)) = Jac_sout_othr_122;%
Jacobian(indices(14),indices(11)) = Jac_sout_othr_222;%
Jacobian(indices(15),indices(12)) = Jac_sout_othr_322;%
Jacobian(indices(16),indices(13)) = Jac_sout_othr_132;%
Jacobian(indices(17),indices(14)) = Jac_sout_othr_232;%
Jacobian(indices(18),indices(15)) = Jac_sout_othr_332;%
Jacobian(indices(19),indices(16)) = Jac_sout_othr_113;%
%Jacobian(indices(20),indices(17)) = Jac_sout_othr_213;%
%Jacobian(indices(21),indices(18)) = Jac_sout_othr_313;%
%Jacobian(indices(22),indices(19)) = Jac_sout_othr_123;%
Jacobian(indices(23),indices(20)) = Jac_sout_othr_223;%
Jacobian(indices(24),indices(21)) = Jac_sout_othr_323;%
Jacobian(indices(25),indices(22)) = Jac_sout_othr_133;%
Jacobian(indices(26),indices(23)) = Jac_sout_othr_233;%
Jacobian(indices(27),indices(24)) = Jac_sout_othr_333;%
%;%
Jacobian(indices(1),indices(10)) = Jac_uppr_othr_111;%
Jacobian(indices(2),indices(11)) = Jac_uppr_othr_211;%
Jacobian(indices(3),indices(12)) = Jac_uppr_othr_311;%
Jacobian(indices(4),indices(13)) = Jac_uppr_othr_121;%
Jacobian(indices(5),indices(14)) = Jac_uppr_othr_221;%
Jacobian(indices(6),indices(15)) = Jac_uppr_othr_321;%
Jacobian(indices(7),indices(16)) = Jac_uppr_othr_131;%
Jacobian(indices(8),indices(17)) = Jac_uppr_othr_231;%
Jacobian(indices(9),indices(18)) = Jac_uppr_othr_331;%
Jacobian(indices(10),indices(19)) = Jac_uppr_othr_112;%
Jacobian(indices(11),indices(20)) = Jac_uppr_othr_212;%
Jacobian(indices(12),indices(21)) = Jac_uppr_othr_312;%
Jacobian(indices(13),indices(22)) = Jac_uppr_othr_122;%
Jacobian(indices(14),indices(23)) = Jac_uppr_othr_222;%
Jacobian(indices(15),indices(24)) = Jac_uppr_othr_322;%
Jacobian(indices(16),indices(25)) = Jac_uppr_othr_132;%
Jacobian(indices(17),indices(26)) = Jac_uppr_othr_232;%
Jacobian(indices(18),indices(27)) = Jac_uppr_othr_332;%
%Jacobian(indices(19),indices(28)) = Jac_uppr_othr_113;%
%Jacobian(indices(20),indices(29)) = Jac_uppr_othr_213;%
%Jacobian(indices(21),indices(30)) = Jac_uppr_othr_313;%
%Jacobian(indices(22),indices(31)) = Jac_uppr_othr_123;%
%Jacobian(indices(23),indices(32)) = Jac_uppr_othr_223;%
%Jacobian(indices(24),indices(33)) = Jac_uppr_othr_323;%
%Jacobian(indices(25),indices(34)) = Jac_uppr_othr_133;%
%Jacobian(indices(26),indices(35)) = Jac_uppr_othr_233;%
%;%
%Jacobian(indices(1),indices(-8)) = Jac_down_othr_111;%
%Jacobian(indices(2),indices(-7)) = Jac_down_othr_211;%
%Jacobian(indices(3),indices(-6)) = Jac_down_othr_311;%
%Jacobian(indices(4),indices(-5)) = Jac_down_othr_121;%
%Jacobian(indices(5),indices(-4)) = Jac_down_othr_221;%
%Jacobian(indices(6),indices(-3)) = Jac_down_othr_321;%
%Jacobian(indices(7),indices(-2)) = Jac_down_othr_131;%
%Jacobian(indices(8),indices(-1)) = Jac_down_othr_231;%
%Jacobian(indices(9),indices(0)) = Jac_down_othr_331;%
Jacobian(indices(10),indices(1)) = Jac_down_othr_112;%
Jacobian(indices(11),indices(2)) = Jac_down_othr_212;%
Jacobian(indices(12),indices(3)) = Jac_down_othr_312;%
Jacobian(indices(13),indices(4)) = Jac_down_othr_122;%
Jacobian(indices(14),indices(5)) = Jac_down_othr_222;%
Jacobian(indices(15),indices(6)) = Jac_down_othr_322;%
Jacobian(indices(16),indices(7)) = Jac_down_othr_132;%
Jacobian(indices(17),indices(8)) = Jac_down_othr_232;%
Jacobian(indices(18),indices(9)) = Jac_down_othr_332;%
Jacobian(indices(19),indices(10)) = Jac_down_othr_113;%
Jacobian(indices(20),indices(11)) = Jac_down_othr_213;%
Jacobian(indices(21),indices(12)) = Jac_down_othr_313;%
Jacobian(indices(22),indices(13)) = Jac_down_othr_123;%
Jacobian(indices(23),indices(14)) = Jac_down_othr_223;%
Jacobian(indices(24),indices(15)) = Jac_down_othr_323;%
Jacobian(indices(25),indices(16)) = Jac_down_othr_133;%
Jacobian(indices(26),indices(17)) = Jac_down_othr_233;%
Jacobian(indices(27),indices(18)) = Jac_down_othr_333;%
                                



end



