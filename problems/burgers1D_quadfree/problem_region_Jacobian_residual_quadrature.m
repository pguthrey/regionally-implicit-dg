function [Jacobian,residual] ...
        = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,data)
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

maxspeedF = 0;

vectpsi_west = data.vectpsi_west;
vectpsi_east = data.vectpsi_east;
rx = 2;
q1 = DGregion(:,rx-1);
q2 = DGregion(:,rx+0);
q3 = DGregion(:,rx+1);

q1_past = DGregion_past(:,rx-1);
q2_past = DGregion_past(:,rx+0);
q3_past = DGregion_past(:,rx+1);

psikpsim_ee_whole = data.psikpsim_ee;
psikpsim_we_whole = data.psikpsim_we;
psikpsim_ew_whole = data.psikpsim_ew;
psikpsim_ww_whole = data.psikpsim_ww;



%RESIDUAL

east_flux_1 = 0;
east_flux_2 = 0;
east_trunc_3 = 0;
west_flux_3 = 0;
west_flux_2 = 0;
west_trunc_1 = 0;

% JACOBIAN 
Jac_east_flux_2 = 0;
Jac_east_flux_1 = 0;
Jac_east_othr_2 = 0;
Jac_east_othr_1 = 0;
Jac_east_trunc_3 = 0;
Jac_west_flux_3 = 0;
Jac_west_flux_2 = 0;
Jac_west_othr_3 = 0;
Jac_west_othr_2 = 0;
Jac_west_trunc_1 = 0;

[Jac_cell_1,residual_cell_1] = problem_jac_res_init(q1,q1_past,data);
[Jac_cell_2,residual_cell_2] = problem_jac_res_init(q2,q2_past,data);
[Jac_cell_3,residual_cell_3] = problem_jac_res_init(q3,q3_past,data);

v1quad = 1;
for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    if space_dims >= 2
        v1quad = Dlist(k,2);
    end    
    weight = Dquadwgts(k);    
    wgtnuv1 = weight*nuv1;
    
    psie = vectpsi_east(:,:,tquad,v1quad);
    psiw = vectpsi_west(:,:,tquad,v1quad);
    
    psieT = psie';
    psiwT = psiw';
    
    psikpsim_ee = psikpsim_ee_whole(:,:,tquad,v1quad);
    psikpsim_we = psikpsim_we_whole(:,:,tquad,v1quad);
    psikpsim_ew = psikpsim_ew_whole(:,:,tquad,v1quad);
    psikpsim_ww = psikpsim_ww_whole(:,:,tquad,v1quad);
                
    % EAST AND WEST TERMS
    w1_p = psie*q1;
    w1_m = psiw*q1;
    
    w2_p = psie*q2;
    w2_m = psiw*q2;
    
    w3_p = psie*q3;
    w3_m = psiw*q3;    
    
    F1_p = w1_p^2/2;
    F2_p = w2_p^2/2;
    F3_p = w3_p^2/2;
    
    F1_m = w1_m^2/2;
    F2_m = w2_m^2/2;
    F3_m = w3_m^2/2;
         
    DF1_pDq = w1_p;
    DF2_pDq = w2_p;
    DF3_pDq = w3_p;
    
    DF1_mDq = w1_m;
    DF2_mDq = w2_m;
    DF3_mDq = w3_m;
    
    
    %psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    east_trunc_3 = east_trunc_3 + psieT*(F3_p*wgtnuv1);
    
    Jac_east_trunc_3 = Jac_east_trunc_3 + psikpsim_ee*(DF3_pDq*wgtnuv1);
    
    west_trunc_1 = west_trunc_1 - psiwT*(F1_m*wgtnuv1);
    
    Jac_west_trunc_1 = Jac_west_trunc_1 - psikpsim_ww*(DF1_mDq*wgtnuv1);
          
    % East terms
    % Interface 1-2
    speed = max(abs(w1_p),abs(w2_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F1_p + F2_m - speed*(w2_m-w1_p))*wgtnuv1/2;
    east_flux_1 = east_flux_1 + psieT*Flux;
    F_dql = (DF1_pDq + speed)*wgtnuv1/2;
    Jac_east_flux_1 = Jac_east_flux_1 + psikpsim_ee*F_dql; 

    F_dqr = (DF2_mDq - speed)*wgtnuv1/2;
    Jac_east_othr_1 = Jac_east_othr_1 + psikpsim_ew*F_dqr;    
    
    % Interface 2-3
    speed = max(abs(w2_p),abs(w3_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F2_p + F3_m - speed*(w3_m-w2_p))*wgtnuv1/2;
    east_flux_2 = east_flux_2 + psieT*Flux;
    F_dql = (DF2_pDq + speed)*wgtnuv1/2;
    Jac_east_flux_2 = Jac_east_flux_2 + psikpsim_ee*F_dql;     
    
    F_dqr = (DF3_mDq - speed)*wgtnuv1/2;
    Jac_east_othr_2 = Jac_east_othr_2 + psikpsim_ew*F_dqr; 
    
    speed = max(abs(w1_p),abs(w2_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F1_p + F2_m - speed*(w2_m-w1_p))*wgtnuv1/2;
    west_flux_2 = west_flux_2 - psiwT*Flux;
    F_dql = (DF1_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF2_mDq - speed)*wgtnuv1/2;

    %psi_m = psiw;
    Jac_west_flux_2 = Jac_west_flux_2 - psikpsim_ww*F_dqr; 

    %psi_m = psie;
    Jac_west_othr_2 = Jac_west_othr_2 - psikpsim_we*F_dql; 
    
    
    % Interface 2-3
    speed = max(abs(w2_p),abs(w3_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F2_p + F3_m - speed*(w3_m-w2_p))*wgtnuv1/2;
    west_flux_3 = west_flux_3 - psiwT*Flux;
    F_dql = (DF2_pDq + speed)*wgtnuv1/2;
    F_dqr = (DF3_mDq - speed)*wgtnuv1/2;
    Jac_west_flux_3 = Jac_west_flux_3 - psikpsim_ww*F_dqr; 

    Jac_west_othr_3 = Jac_west_othr_3 - psikpsim_we*F_dql; 
   
    
    
end

% RESIDUAL INFORMATION
indices = @(i) thetaT*(i-1) + (1:thetaT);

Jacobian(indices(1),indices(1)) =   Jac_cell_1 + ...
                                    Jac_east_flux_1 + ... 
                                    Jac_west_trunc_1;
Jacobian(indices(1),indices(2)) = Jac_east_othr_1;


Jacobian(indices(2),indices(1)) = Jac_west_othr_2;
Jacobian(indices(2),indices(2)) =   Jac_cell_2 + ...
                                    Jac_east_flux_2 + ... 
                                    Jac_west_flux_2;
Jacobian(indices(2),indices(3)) = Jac_east_othr_2;


Jacobian(indices(3),indices(2)) = Jac_west_othr_3;
Jacobian(indices(3),indices(3)) =   Jac_cell_3 + ...
                                    Jac_east_trunc_3 + ... 
                                    Jac_west_flux_3;
                                
                                         
                                
residual(indices(1),1) =   residual_cell_1 + ...
                                    east_flux_1 + ... 
                                    west_trunc_1;
                                
residual(indices(2),1) =   residual_cell_2 + ...
                                    east_flux_2 + ... 
                                    west_flux_2 ;
residual(indices(3),1) =    residual_cell_3 + ...
                                    east_trunc_3 + ... 
                                    west_flux_3 ;
                                
    
end



