function [east_trunc,west_trunc, ...
          east_flux, west_flux, ...
          Jac_east_cell,Jac_west_cell, ...
        Jac_east_othr,Jac_west_othr, ...
        Jac_east_trunc,Jac_west_trunc, ...
        nort_trunc,sout_trunc, ...
          nort_flux, sout_flux, ...
          Jac_nort_cell,Jac_sout_cell, ...
        Jac_nort_othr,Jac_sout_othr, ...
        Jac_nort_trunc,Jac_sout_trunc, ...
        maxspeedF,maxspeedG] ...
        = problem_precompute_region_all(DGregion,data,~)
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


Jac_east_cell = zeros(thetaT,thetaT,3,3);
Jac_west_cell = zeros(thetaT,thetaT,3,3);    
Jac_east_othr = zeros(thetaT,thetaT,3,3);
Jac_west_othr = zeros(thetaT,thetaT,3,3);  
Jac_east_trunc = zeros(thetaT,thetaT,3,3);   
Jac_west_trunc = zeros(thetaT,thetaT,3,3); 

east_trunc = zeros(thetaT,3,3); 
west_trunc = zeros(thetaT,3,3); 
east_flux = zeros(thetaT,3,3); 
west_flux = zeros(thetaT,3,3); 

nort_trunc = zeros(thetaT,3,3); 
sout_trunc = zeros(thetaT,3,3); 
nort_flux = zeros(thetaT,3,3); 
sout_flux = zeros(thetaT,3,3); 

Jac_nort_cell = zeros(thetaT,thetaT,3,3);
Jac_sout_cell = zeros(thetaT,thetaT,3,3);    
Jac_nort_othr = zeros(thetaT,thetaT,3,3);
Jac_sout_othr = zeros(thetaT,thetaT,3,3);  
Jac_nort_trunc = zeros(thetaT,thetaT,3,3);   
Jac_sout_trunc = zeros(thetaT,thetaT,3,3); 


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


%{
v2quad = 1;
v3quad = 1;
Jac_cell = zeros(thetaT,thetaT,3,3);
tic
for rx = 1:3
    for ry = 1:3    
        qstar = DGregion(:,rx,ry);
        Jac_cell_x = problem_exact_jaccell_eastwest_M4(qstar,nuv1); 
        Jac_cell_y = problem_exact_jaccell_nortsout_M4(qstar,nuv2);        
        Jac_cell(:,:,rx,ry) =  Jac_cell_x + Jac_cell_y;
    end
end
disp(['symbolic 1 = ' num2str(toc)])
%}



%disp(['symbolic 2 = ' num2str(toc)])

%{
%Jac_cell = zeros(thetaT,thetaT,3,3);
%tic
for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        v1quad = Dp1list(k,2); 
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            if space_dims >= 3
                v3quad = Dp1list(k,4);
            end
        end
        weight = Dp1quadwgts(k);

        psidxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        psideta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        psi = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        
        %psidxiipsi = psidxii*psi;
        %psidetapsi = psideta*psi;
        
        for rx = 1:3
            for ry = 1:3
                %Derivative terms
                qstar = DGregion(:,rx,ry);
                wstar = psi*qstar;
                Fstar = wstar^2/2;
                fterm = nuv1*weight*psidxii*Fstar;
                gterm = nuv2*weight*psideta*Fstar;        
                                
                residual_cell(:,rx,ry) = residual_cell(:,rx,ry) - fterm - gterm; 
                maxspeedF = max(maxspeedF,abs(wstar));
                maxspeedG = max(maxspeedG,abs(wstar));

                %psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
                %psim = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
                %fterm = weight*psidxiipsi*wstar;
                %gterm = weight*psidetapsi*wstar;
                %Jac_cell(:,:,rx,ry) = Jac_cell(:,:,rx,ry) - fterm - gterm;
            end
        end
end 
%disp(['brute 1 = ' num2str(toc)])

%}
%{
residual_cell = zeros(thetaT,3,3);
Jac_cell = zeros(thetaT,thetaT,3,3);
tic 
for rx = 1:3
    for ry = 1:3
        qstar = DGregion(:,rx,ry);
        residual_cell(:,rx,ry) = tauflux*qstar-Ipast*DGregion_past(:,rx,ry);
        Jac_cell(:,:,rx,ry) = tauflux;
        res_cell_this = zeros(thetaT,1);
        Jac_cell_this = zeros(thetaT,thetaT);
        for k = 1:Pdp1
            tquad = Dp1list(k,1); 
            v1quad = Dp1list(k,2); 
            if space_dims >= 2
                v2quad = Dp1list(k,3);
                if space_dims >= 3
                    v3quad = Dp1list(k,4);
                end
            end
            weight = Dp1quadwgts(k);

            psidxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
            psideta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
            psi = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);

            psidxiipsi = psidxii*psi;
            psidetapsi = psideta*psi;

            %Derivative terms
            wstar = psi*qstar;
            Fstar = wstar^2/2;
            fterm = nuv1*weight*psidxii*Fstar;
            gterm = nuv2*weight*psideta*Fstar;        
            res_cell_this = res_cell_this - fterm;
            res_cell_this = res_cell_this - gterm; 
            maxspeedF = max(maxspeedF,abs(wstar));
            maxspeedG = max(maxspeedG,abs(wstar));

            %psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
            %psim = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
            fterm = nuv1*weight*psidxiipsi*wstar;
            gterm = nuv2*weight*psidetapsi*wstar;
            Jac_cell_this = Jac_cell_this - fterm;
            Jac_cell_this = Jac_cell_this - gterm;
        end
        residual_cell(:,rx,ry) = residual_cell(:,rx,ry) + res_cell_this; 
        Jac_cell(:,:,rx,ry) = Jac_cell(:,:,rx,ry) + Jac_cell_this;
    end
end
disp(['brute 2 = ' num2str(toc)])
disp('-----------------------------------')
%}

v2quad = 1;
%v3quad = 1;

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


psikpsim_ee_whole = data.psikpsim_ee;
psikpsim_we_whole = data.psikpsim_we;
psikpsim_ew_whole = data.psikpsim_ew;
psikpsim_ww_whole = data.psikpsim_ww;

psikpsim_nn_whole = data.psikpsim_nn;
psikpsim_sn_whole = data.psikpsim_sn;
psikpsim_ns_whole = data.psikpsim_ns;
psikpsim_ss_whole = data.psikpsim_ss;


Jac_sout_othr_33 = zeros(thetaT,thetaT);


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
    
    psie = vectpsi_east(:,:,tquad,v1quad,v2quad);
    psiw = vectpsi_west(:,:,tquad,v1quad,v2quad);
    psin = vectpsi_nort(:,:,tquad,v1quad,v2quad);
    psis = vectpsi_sout(:,:,tquad,v1quad,v2quad);
    
    psieT = psie';
    psiwT = psiw';
    psinT = psin';
    psisT = psis';
    
    psikpsim_ee = psikpsim_ee_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_we = psikpsim_we_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_ew = psikpsim_ew_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_ww = psikpsim_ww_whole(:,:,tquad,v1quad,v2quad);
    
    psikpsim_nn = psikpsim_nn_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_sn = psikpsim_sn_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_ns = psikpsim_ns_whole(:,:,tquad,v1quad,v2quad);
    psikpsim_ss = psikpsim_ss_whole(:,:,tquad,v1quad,v2quad);   
    
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
    
    %psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    east_trunc(:,1,1) = east_trunc(:,1,1) + (psieT*F1_p)*(nuv1*weight);
    east_trunc(:,2,1) = east_trunc(:,2,1) + (psieT*F2_p)*(nuv1*weight);
    east_trunc(:,3,1) = east_trunc(:,3,1) + (psieT*F3_p)*(nuv1*weight);
    east_trunc(:,1,2) = east_trunc(:,1,2) + (psieT*F4_p)*(nuv1*weight);
    east_trunc(:,2,2) = east_trunc(:,2,2) + (psieT*F5_p)*(nuv1*weight);
    east_trunc(:,3,2) = east_trunc(:,3,2) + (psieT*F6_p)*(nuv1*weight);
    east_trunc(:,1,3) = east_trunc(:,1,3) + (psieT*F7_p)*(nuv1*weight);
    east_trunc(:,2,3) = east_trunc(:,2,3) + (psieT*F8_p)*(nuv1*weight);
    east_trunc(:,3,3) = east_trunc(:,3,3) + (psieT*F9_p)*(nuv1*weight);

    %psim = psie;
    %psikpsim = psik*psim;    
    Jac_east_trunc_11 = Jac_east_trunc_11 + (psikpsim_ee*DF1_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,2,1) = Jac_east_trunc(:,:,2,1) + (psikpsim_ee*DF2_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,3,1) = Jac_east_trunc(:,:,3,1) + (psikpsim_ee*DF3_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,1,2) = Jac_east_trunc(:,:,1,2) + (psikpsim_ee*DF4_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,2,2) = Jac_east_trunc(:,:,2,2) + (psikpsim_ee*DF5_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,3,2) = Jac_east_trunc(:,:,3,2) + (psikpsim_ee*DF6_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,1,3) = Jac_east_trunc(:,:,1,3) + (psikpsim_ee*DF7_pDq)*(nuv1*weight);
    Jac_east_trunc(:,:,2,3) = Jac_east_trunc(:,:,2,3) + (psikpsim_ee*DF8_pDq)*(nuv1*weight);
    Jac_east_trunc_33 = Jac_east_trunc_33 + (psikpsim_ee*DF9_pDq)*(nuv1*weight);      
    
    %psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    west_trunc(:,1,1) = west_trunc(:,1,1) - (psiwT*F1_m)*(nuv1*weight);
    west_trunc(:,2,1) = west_trunc(:,2,1) - (psiwT*F2_m)*(nuv1*weight);
    west_trunc(:,3,1) = west_trunc(:,3,1) - (psiwT*F3_m)*(nuv1*weight);
    west_trunc(:,1,2) = west_trunc(:,1,2) - (psiwT*F4_m)*(nuv1*weight);
    west_trunc(:,2,2) = west_trunc(:,2,2) - (psiwT*F5_m)*(nuv1*weight);
    west_trunc(:,3,2) = west_trunc(:,3,2) - (psiwT*F6_m)*(nuv1*weight);
    west_trunc(:,1,3) = west_trunc(:,1,3) - (psiwT*F7_m)*(nuv1*weight);
    west_trunc(:,2,3) = west_trunc(:,2,3) - (psiwT*F8_m)*(nuv1*weight);
    west_trunc(:,3,3) = west_trunc(:,3,3) - (psiwT*F9_m)*(nuv1*weight);  
    
    %psim = psiw;
    %psikpsim = psik*psim;    
    Jac_west_trunc_11 = Jac_west_trunc_11 - (psikpsim_ww*DF1_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,2,1) = Jac_west_trunc(:,:,2,1) - (psikpsim_ww*DF2_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,3,1) = Jac_west_trunc(:,:,3,1) - (psikpsim_ww*DF3_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,1,2) = Jac_west_trunc(:,:,1,2) - (psikpsim_ww*DF4_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,2,2) = Jac_west_trunc(:,:,2,2) - (psikpsim_ww*DF5_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,3,2) = Jac_west_trunc(:,:,3,2) - (psikpsim_ww*DF6_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,1,3) = Jac_west_trunc(:,:,1,3) - (psikpsim_ww*DF7_mDq)*(nuv1*weight);
    Jac_west_trunc(:,:,2,3) = Jac_west_trunc(:,:,2,3) - (psikpsim_ww*DF8_mDq)*(nuv1*weight);
    Jac_west_trunc_33 = Jac_west_trunc_33 - (psikpsim_ww*DF9_mDq)*(nuv1*weight);     
          
    % East terms
    %psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    % Interface 1-2
    speed = max(abs(w1_p),abs(w2_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F1_p + F2_m - speed*(w2_m-w1_p))/2;
    east_flux(:,1,1) = east_flux(:,1,1) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF1_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell_11 = Jac_east_cell_11 + dFlux;    

    F_dqr = (DF2_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ew*F_dqr.*weight*nuv1;
    Jac_east_othr_11 = Jac_east_othr_11 + dFlux;    
    
    
    % Interface 2-3
    speed = max(abs(w2_p),abs(w3_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F2_p + F3_m - speed*(w3_m-w2_p))/2;
    east_flux(:,2,1) = east_flux(:,2,1) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF2_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell(:,:,2,1) = Jac_east_cell(:,:,2,1) + dFlux;    

    F_dqr = (DF3_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ew*F_dqr.*weight*nuv1;
    Jac_east_othr(:,:,2,1) = Jac_east_othr(:,:,2,1) + dFlux;    

    
    % Interface 4-5
    speed = max(abs(w4_p),abs(w5_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F4_p + F5_m - speed*(w5_m-w4_p))/2;
    east_flux(:,1,2) = east_flux(:,1,2) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF4_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell(:,:,1,2) = Jac_east_cell(:,:,1,2) + dFlux;    

    F_dqr = (DF5_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ew*F_dqr.*weight*nuv1;
    Jac_east_othr(:,:,1,2) = Jac_east_othr(:,:,1,2) + dFlux;    
    
    
    % Interface 5-6
    speed = max(abs(w5_p),abs(w6_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F5_p + F6_m - speed*(w6_m-w5_p))/2;
    east_flux(:,2,2) = east_flux(:,2,2) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF5_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell(:,:,2,2) = Jac_east_cell(:,:,2,2) + dFlux;    

    F_dqr = (DF6_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ew*F_dqr.*weight*nuv1;
    Jac_east_othr(:,:,2,2) = Jac_east_othr(:,:,2,2) + dFlux;    

    % Interface 7-8
    speed = max(abs(w7_p),abs(w8_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F7_p + F8_m - speed*(w8_m-w7_p))/2;
    east_flux(:,1,3) = east_flux(:,1,3) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF7_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell(:,:,1,3) = Jac_east_cell(:,:,1,3) + dFlux;    

    F_dqr = (DF8_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ew*F_dqr.*weight*nuv1;
    Jac_east_othr(:,:,1,3) = Jac_east_othr(:,:,1,3) + dFlux;    
    
    % Interface 8-9
    speed = max(abs(w8_p),abs(w9_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F8_p + F9_m - speed*(w9_m-w8_p))/2;
    east_flux(:,2,3) = east_flux(:,2,3) + (psieT*Flux)*(nuv1*weight);
    F_dql = (DF8_pDq + speed)/2;
    %psi_m = psie;
    dFlux = psikpsim_ee*F_dql.*weight*nuv1;
    Jac_east_cell(:,:,2,3) = Jac_east_cell(:,:,2,3) + dFlux;    

    F_dqr = (DF9_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ee*F_dqr.*weight*nuv1;
    Jac_east_othr(:,:,2,3) = Jac_east_othr(:,:,2,3) + dFlux;    
    
    % West terms
    %psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    % Interface 1-2
    speed = max(abs(w1_p),abs(w2_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F1_p + F2_m - speed*(w2_m-w1_p))/2;
    west_flux(:,2,1) = west_flux(:,2,1) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF1_pDq + speed)/2;
    F_dqr = (DF2_mDq - speed)/2;

    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell(:,:,2,1) = Jac_west_cell(:,:,2,1) - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr(:,:,2,1) = Jac_west_othr(:,:,2,1) - dFlux;    
    
    
    % Interface 2-3
    speed = max(abs(w2_p),abs(w3_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F2_p + F3_m - speed*(w3_m-w2_p))/2;
    west_flux(:,3,1) = west_flux(:,3,1) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF2_pDq + speed)/2;
    F_dqr = (DF3_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell(:,:,3,1) = Jac_west_cell(:,:,3,1) - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr(:,:,3,1) = Jac_west_othr(:,:,3,1) - dFlux;    

    
    % Interface 4-5
    speed = max(abs(w4_p),abs(w5_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F4_p + F5_m - speed*(w5_m-w4_p))/2;
    west_flux(:,2,2) = west_flux(:,2,2) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF4_pDq + speed)/2;
    F_dqr = (DF5_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell(:,:,2,2) = Jac_west_cell(:,:,2,2) - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr(:,:,2,2) = Jac_west_othr(:,:,2,2) - dFlux;    
        
    % Interface 5-6
    speed = max(abs(w5_p),abs(w6_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F5_p + F6_m - speed*(w6_m-w5_p))/2;
    west_flux(:,3,2) = west_flux(:,3,2) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF5_pDq + speed)/2;
    F_dqr = (DF6_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell(:,:,3,2) = Jac_west_cell(:,:,3,2) - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr(:,:,3,2) = Jac_west_othr(:,:,3,2) - dFlux;    

    % Interface 7-8
    speed = max(abs(w7_p),abs(w8_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F7_p + F8_m - speed*(w8_m-w7_p))/2;
    west_flux(:,2,3) = west_flux(:,2,3) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF7_pDq + speed)/2;
    F_dqr = (DF8_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell(:,:,2,3) = Jac_west_cell(:,:,2,3) - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr(:,:,2,3) = Jac_west_othr(:,:,2,3) - dFlux;    
    
    % Interface 8-9
    speed = max(abs(w8_p),abs(w9_m));
    maxspeedF = max(maxspeedF,speed);
    Flux = (F8_p + F9_m - speed*(w9_m-w8_p))/2;
    west_flux(:,3,3) = west_flux(:,3,3) - (psiwT*Flux)*(nuv1*weight);
    F_dql = (DF8_pDq + speed)/2;
    F_dqr = (DF9_mDq - speed)/2;
    %psi_m = psiw;
    dFlux = psikpsim_ww*F_dqr.*weight*nuv1;
    Jac_west_cell_33 = Jac_west_cell_33 - dFlux;    

    %psi_m = psie;
    dFlux = psikpsim_we*F_dql.*weight*nuv1;
    Jac_west_othr_33 = Jac_west_othr_33 - dFlux;    
    
    

    %NORT SOUT terms --------------------------------------------\
    %=============================================================
    w1_p = psin*q1;
    w1_m = psis*q1;

    w2_p = psin*q2;
    w2_m = psis*q2;    

    w3_p = psin*q3;
    w3_m = psis*q3;    

    w4_p = psin*q4;
    w4_m = psis*q4;    

    w5_p = psin*q5;
    w5_m = psis*q5;    

    w6_p = psin*q6;
    w6_m = psis*q6;    

    w7_p = psin*q7;
    w7_m = psis*q7;    

    w8_p = psin*q8;
    w8_m = psis*q8;    

    w9_p = psin*q9;
    w9_m = psis*q9;    
    
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
    
    %%{
    %psik = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    nort_trunc(:,1,1) = nort_trunc(:,1,1) + (psinT*F1_p)*(nuv2*weight);
    nort_trunc(:,2,1) = nort_trunc(:,2,1) + (psinT*F2_p)*(nuv2*weight);
    nort_trunc(:,3,1) = nort_trunc(:,3,1) + (psinT*F3_p)*(nuv2*weight);
    nort_trunc(:,1,2) = nort_trunc(:,1,2) + (psinT*F4_p)*(nuv2*weight);
    nort_trunc(:,2,2) = nort_trunc(:,2,2) + (psinT*F5_p)*(nuv2*weight);
    nort_trunc(:,3,2) = nort_trunc(:,3,2) + (psinT*F6_p)*(nuv2*weight);
    nort_trunc(:,1,3) = nort_trunc(:,1,3) + (psinT*F7_p)*(nuv2*weight);
    nort_trunc(:,2,3) = nort_trunc(:,2,3) + (psinT*F8_p)*(nuv2*weight);
    nort_trunc(:,3,3) = nort_trunc(:,3,3) + (psinT*F9_p)*(nuv2*weight);

    %psim = psin;
    %psikpsim = psik*psim;
    Jac_nort_trunc_11 = Jac_nort_trunc_11 + (psikpsim_nn*DF1_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,2,1) = Jac_nort_trunc(:,:,2,1) + (psikpsim_nn*DF2_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,3,1) = Jac_nort_trunc(:,:,3,1) + (psikpsim_nn*DF3_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,1,2) = Jac_nort_trunc(:,:,1,2) + (psikpsim_nn*DF4_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,2,2) = Jac_nort_trunc(:,:,2,2) + (psikpsim_nn*DF5_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,3,2) = Jac_nort_trunc(:,:,3,2) + (psikpsim_nn*DF6_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,1,3) = Jac_nort_trunc(:,:,1,3) + (psikpsim_nn*DF7_pDq)*(nuv2*weight);
    Jac_nort_trunc(:,:,2,3) = Jac_nort_trunc(:,:,2,3) + (psikpsim_nn*DF8_pDq)*(nuv2*weight);
    Jac_nort_trunc_33 = Jac_nort_trunc_33 + (psikpsim_nn*DF9_pDq)*(nuv2*weight);      
    
    %psik = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    sout_trunc(:,1,1) = sout_trunc(:,1,1) - (psisT*F1_m)*(nuv2*weight);
    sout_trunc(:,2,1) = sout_trunc(:,2,1) - (psisT*F2_m)*(nuv2*weight);
    sout_trunc(:,3,1) = sout_trunc(:,3,1) - (psisT*F3_m)*(nuv2*weight);
    sout_trunc(:,1,2) = sout_trunc(:,1,2) - (psisT*F4_m)*(nuv2*weight);
    sout_trunc(:,2,2) = sout_trunc(:,2,2) - (psisT*F5_m)*(nuv2*weight);
    sout_trunc(:,3,2) = sout_trunc(:,3,2) - (psisT*F6_m)*(nuv2*weight);
    sout_trunc(:,1,3) = sout_trunc(:,1,3) - (psisT*F7_m)*(nuv2*weight);
    sout_trunc(:,2,3) = sout_trunc(:,2,3) - (psisT*F8_m)*(nuv2*weight);
    sout_trunc(:,3,3) = sout_trunc(:,3,3) - (psisT*F9_m)*(nuv2*weight);  
    
    %psim = psis;
    %psikpsim = psik*psim;
    Jac_sout_trunc_11 = Jac_sout_trunc_11 - (psikpsim_ss*DF1_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,2,1) = Jac_sout_trunc(:,:,2,1) - (psikpsim_ss*DF2_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,3,1) = Jac_sout_trunc(:,:,3,1) - (psikpsim_ss*DF3_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,1,2) = Jac_sout_trunc(:,:,1,2) - (psikpsim_ss*DF4_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,2,2) = Jac_sout_trunc(:,:,2,2) - (psikpsim_ss*DF5_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,3,2) = Jac_sout_trunc(:,:,3,2) - (psikpsim_ss*DF6_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,1,3) = Jac_sout_trunc(:,:,1,3) - (psikpsim_ss*DF7_mDq)*(nuv2*weight);
    Jac_sout_trunc(:,:,2,3) = Jac_sout_trunc(:,:,2,3) - (psikpsim_ss*DF8_mDq)*(nuv2*weight);
    Jac_sout_trunc_33 = Jac_sout_trunc_33 - (psikpsim_ss*DF9_mDq)*(nuv2*weight);           
    %}
    
    % nort terms -------------------------------------------------------
    %psik = vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
    % Interface 1-4
    speed = max(abs(w1_p),abs(w4_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F1_p + F4_m - speed*(w4_m-w1_p))/2;
    nort_flux(:,1,1) = nort_flux(:,1,1) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF1_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell_11 = Jac_nort_cell_11 + dFlux;    

    F_dqr = (DF4_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr_11 = Jac_nort_othr_11 + dFlux;    
    
    % Interface 2-5
    speed = max(abs(w2_p),abs(w5_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F2_p + F5_m - speed*(w5_m-w2_p))/2;
    nort_flux(:,2,1) = nort_flux(:,2,1) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF2_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell(:,:,2,1) = Jac_nort_cell(:,:,2,1) + dFlux;    

    F_dqr = (DF5_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr(:,:,2,1) = Jac_nort_othr(:,:,2,1) + dFlux;  
    
    % Interface 3-6
    speed = max(abs(w3_p),abs(w6_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F3_p + F6_m - speed*(w6_m-w3_p))/2;
    nort_flux(:,3,1) = nort_flux(:,3,1) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF3_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell(:,:,3,1) = Jac_nort_cell(:,:,3,1) + dFlux;    

    F_dqr = (DF6_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr(:,:,3,1) = Jac_nort_othr(:,:,3,1) + dFlux;     
    
    % Interface 4-7
    speed = max(abs(w4_p),abs(w7_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F4_p + F7_m - speed*(w7_m-w4_p))/2;
    nort_flux(:,1,2) = nort_flux(:,1,2) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF4_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell(:,:,1,2) = Jac_nort_cell(:,:,1,2) + dFlux;    

    F_dqr = (DF7_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr(:,:,1,2) = Jac_nort_othr(:,:,1,2) + dFlux;    
    
    % Interface 5-8
    speed = max(abs(w5_p),abs(w8_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F5_p + F8_m - speed*(w8_m-w5_p))/2;
    nort_flux(:,2,2) = nort_flux(:,2,2) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF5_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell(:,:,2,2) = Jac_nort_cell(:,:,2,2) + dFlux;    

    F_dqr = (DF8_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr(:,:,2,2) = Jac_nort_othr(:,:,2,2) + dFlux;  
    
    % Interface 6-9
    speed = max(abs(w6_p),abs(w9_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F6_p + F9_m - speed*(w9_m-w6_p))/2;
    nort_flux(:,3,2) = nort_flux(:,3,2) + (psinT*Flux)*(nuv2*weight);
    F_dql = (DF6_pDq + speed)/2;
    %psi_m = psin;
    dFlux = psikpsim_nn*F_dql.*weight*nuv2;
    Jac_nort_cell(:,:,3,2) = Jac_nort_cell(:,:,3,2) + dFlux;    

    F_dqr = (DF9_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ns*F_dqr.*weight*nuv2;
    Jac_nort_othr(:,:,3,2) = Jac_nort_othr(:,:,3,2) + dFlux;     
    
    % SOUT TERMS -----------------------------------------------------
    % sout terms
    %psik = vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
    % Interface 1-4
    speed = max(abs(w1_p),abs(w4_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F1_p + F4_m - speed*(w4_m-w1_p))/2;
    F_dql = (DF1_pDq + speed)/2;
    F_dqr = (DF4_mDq - speed)/2;

    sout_flux(:,1,2) = sout_flux(:,1,2) - (psisT*Flux)*(nuv2*weight);
    
    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell(:,:,1,2) = Jac_sout_cell(:,:,1,2) - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr(:,:,1,2) = Jac_sout_othr(:,:,1,2) - dFlux;    
    
    % Interface 2-5
    speed = max(abs(w2_p),abs(w5_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F2_p + F5_m - speed*(w5_m-w2_p))/2;
    sout_flux(:,2,2) = sout_flux(:,2,2) - (psisT*Flux)*(nuv2*weight);
    F_dql = (DF2_pDq + speed)/2;
    F_dqr = (DF5_mDq - speed)/2;

    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell(:,:,2,2) = Jac_sout_cell(:,:,2,2) - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr(:,:,2,2) = Jac_sout_othr(:,:,2,2) - dFlux;  
    
    % Interface 3-6
    speed = max(abs(w3_p),abs(w6_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F3_p + F6_m - speed*(w6_m-w3_p))/2;
    sout_flux(:,3,2) = sout_flux(:,3,2) - (psisT*Flux)*(nuv2*weight);
    F_dql = (DF3_pDq + speed)/2;
    F_dqr = (DF6_mDq - speed)/2;

    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell(:,:,3,2) = Jac_sout_cell(:,:,3,2) - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr(:,:,3,2) = Jac_sout_othr(:,:,3,2) - dFlux;     
    
    % Interface 4-7
    speed = max(abs(w4_p),abs(w7_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F4_p + F7_m - speed*(w7_m-w4_p))/2;
    sout_flux(:,1,3) = sout_flux(:,1,3) - (psisT*Flux)*(nuv2*weight);
    F_dql = (DF4_pDq + speed)/2;
    F_dqr = (DF7_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell(:,:,1,3) = Jac_sout_cell(:,:,1,3) - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr(:,:,1,3) = Jac_sout_othr(:,:,1,3) - dFlux;    
    
    % Interface 5-8
    speed = max(abs(w5_p),abs(w8_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F5_p + F8_m - speed*(w8_m-w5_p))/2;
    sout_flux(:,2,3) = sout_flux(:,2,3) - (psisT*Flux)*(nuv2*weight);
    F_dql = (DF5_pDq + speed)/2;
    F_dqr = (DF8_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell(:,:,2,3) = Jac_sout_cell(:,:,2,3) - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr(:,:,2,3) = Jac_sout_othr(:,:,2,3) - dFlux;  
    
    % Interface 6-9
    speed = max(abs(w6_p),abs(w9_m));
    maxspeedG = max(maxspeedG,speed);
    Flux = (F6_p + F9_m - speed*(w9_m-w6_p))/2;
    sout_flux(:,3,3) = sout_flux(:,3,3) - (psisT*Flux)*(nuv2*weight);
    F_dql = (DF6_pDq + speed)/2;
    F_dqr = (DF9_mDq - speed)/2;
    %psi_m = psis;
    dFlux = psikpsim_ss*F_dqr.*weight*nuv2;
    Jac_sout_cell_33 = Jac_sout_cell_33 - dFlux;    

    %psi_m = psin;
    dFlux = psikpsim_sn*F_dql.*weight*nuv2;
    Jac_sout_othr_33 = Jac_sout_othr_33 - dFlux;             
    Jac_sout_othr_33 = Jac_sout_othr_33 - dFlux;         
    
end

Jac_sout_othr_33 = Jac_sout_othr_33;

%%{


end



