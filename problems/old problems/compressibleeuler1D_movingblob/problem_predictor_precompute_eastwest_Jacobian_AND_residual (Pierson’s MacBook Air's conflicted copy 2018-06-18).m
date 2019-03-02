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

thetaT = data.thetaT;

Jac_cell = data.Psi_futr;%data.I_futr - data.Psi_dtau;

Jac_east_cell_temp = zeros(thetaT,thetaT);
Jac_west_cell_temp = zeros(thetaT,thetaT);    

Jac_east_other_temp = zeros(thetaT,thetaT);
Jac_west_other_temp = zeros(thetaT,thetaT);  

Jac_trunc_east_temp = zeros(thetaT,thetaT);   
Jac_trunc_west_temp = zeros(thetaT,thetaT); 

trunceast_temp = 0; 
truncwest_temp = 0; 

fluxeast_temp = 0;
fluxwest_temp = 0;    

residual_cell = data.Psi_futr*qstar-data.I_past*qpast;
thetaTtilde = thetaT/data.Neqns;

%NewtonJacobian_cc = data.NewtonJacobian_cc;
%NewtonJacobian_cc = NewtonJacobian_cc + NewtonJacobian_cc;
%shortcut = NewtonJacobian_cc(1:thetaTtilde,1:thetaTtilde)*qeast(1:thetaTtilde,1);

gmma = data.appdata.gamma;
gm1 = gmma-1;
gm3 = gm1-2;

vectpsi_east = data.vectpsi_east;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west = data.vectpsi_west;
vectpsi_west_Trans = data.vectpsi_west_Trans;
vectpsi = data.vectpsi;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;

nuv1 = data.nuv1;

maxspeedF = 0;
speedeast2 = 0;
speedwest2 = 0;

%{
uvect = 0;
for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        v1quad = data.Dp1list(k,2); 
        weight = data.Dp1quadwgts(k);
        %Derivative terms
        %wstar = vectpsi(1,1:thetaTtilde,tquad,v1quad)*qstar;
        wstar = vectpsi(:,:,tquad,v1quad)*qstar;
        uvect = uvect + vectpsi(1,1:thetaTtilde,tquad,v1quad)'*wstar(2)/wstar(1);        
end
%}
Neqns = data.Neqns;
Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
Dp1list = data.Dp1list;
Dquadwgts = data.Dquadwgts;
Dp1quadwgts = data.Dp1quadwgts;
r_param = data.r_param;

%vectpsi_east_one  = vectpsi_east(1,1:thetaTtilde,:,:);
%rhowest(:,1) = qwest(0*thetaTtilde + (1:thetaTtilde));
%mwest(:,1) = qwest(1*thetaTtilde + (1:thetaTtilde));
%Ewest(:,1) = qwest(2*thetaTtilde + (1:thetaTtilde));

I_mat = eye(Neqns);
for k = 1:Pd
    tquad = Dlist(k,1); 
    weight = Dquadwgts(k);    
    
    wwest = vectpsi_east(:,:,tquad)*qwest;       
    weast = vectpsi_west(:,:,tquad)*qeast;

    %East Boundary terms
    wstar = vectpsi_east(:,:,tquad)*qstar;
       
    rho = wstar(1);
    m = wstar(2);
    E = wstar(3);
    
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    u3 = u2*u;
    p = gm1*(E - rho*u2*0.5);
    h = u*(E+p);
    c = sqrt(gmma*p*orho);
    
    JF_star(3,3) = gmma*u ;
    JF_star(3,2) = gmma*E*orho - gm1*1.5*u2   ;
    JF_star(3,1) = -gmma*u*E*orho + gm1*u3 ;
    JF_star(2,3) = gm1;
    JF_star(2,2) = -gm3*u;
    JF_star(2,1) = gm3*u2*0.5;
    JF_star(1,2) = 1;

    F_star(3,1) = h ;
    F_star(2,1) = m*u + p  ;
    F_star(1,1) = m ;
    
    speedF_star = abs(u)+c;

    dspeeddq_star = 0;        
        %[abs(u) sign(u) 0]*orho + dcdq;    
            
    trunceast_temp = trunceast_temp + (nuv1*weight)*(vectpsi_east_Trans(:,:,tquad)*F_star);
    if r_param > 0
        
        rho = weast(1);
        m = weast(2);
        E = weast(3);

        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        u3 = u2*u;
        p = gm1*(E - rho*u2*0.5);
        h = u*(E+p);
        c = sqrt(gmma*p*orho);

        JF_east(3,3) = gmma*u ;
        JF_east(3,2) = gmma*E*orho - gm1*1.5*u2   ;
        JF_east(3,1) = -gmma*u*E*orho + gm1*u3 ;
        JF_east(2,3) = gm1;
        JF_east(2,2) = -gm3*u;
        JF_east(2,1) = gm3*u2*0.5;
        JF_east(1,2) = 1;

        F_east(3,1) = h ;
        F_east(2,1) = m*u + p  ;
        F_east(1,1) = m ;

        speedF_east = abs(u)+c;

        dspeeddq_east = 0;        
        
        wavg = (weast+wstar)*0.5;
        rho = wavg(1);
        m = wavg(2);
        E = wavg(3);
        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        p = gm1*(E - rho*u2*0.5);        
        c = sqrt(gmma*p*orho);        
        speedF_avg = abs(u)+c;
        
        speedeast = max(speedF_star,speedF_east);
        speedeast = max(speedeast,speedF_avg);
        Flux = (F_star+F_east-speedeast.*(weast-wstar))*0.5;
        fluxeast_temp = fluxeast_temp + (nuv1*weight)*(vectpsi_east_Trans(:,:,tquad)*Flux);  
    end
    
    maxspeedF = max(speedeast2,maxspeedF);

    [ Flux ] = vectpsi_east_Trans(:,:,tquad)*JF_star*vectpsi_east(:,:,tquad).*(weight*nuv1);
    Jac_trunc_east_temp = Jac_trunc_east_temp + Flux;

    if r_param > 0
        dspeeddql = 0;%compute_F_dwavespeed_dql(wstar,weast,quadpoint,1);
        F_dql = (JF_star+speedeast.*I_mat - (weast-wstar)*dspeeddql)*0.5;
        dFlux = vectpsi_east_Trans(:,:,tquad)*F_dql*vectpsi_east(:,:,tquad).*(weight*nuv1);
        Jac_east_cell_temp = Jac_east_cell_temp + dFlux;
        %psi_m = vectpsi_west(:,:,tquad);
        dspeeddqr = 0;%compute_F_dwavespeed_dqr(wstar,weast,quadpoint,1);        
        F_dqr = (JF_east-speedeast.*I_mat - (weast-wstar)*dspeeddqr)*0.5; 
        %????
        dFlux = vectpsi_east_Trans(:,:,tquad)*F_dqr*vectpsi_west(:,:,tquad).*(weight*nuv1);
        Jac_east_other_temp = Jac_east_other_temp + dFlux;
    end
    %west boundary terms
    wstar = vectpsi_west(:,:,tquad)*qstar;
    rho = wstar(1);
    m = wstar(2);
    E = wstar(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    speedwest1 = abs(u) + c;

    wstar = vectpsi_west(:,:,tquad)*qstar;
    
    
    rho = wstar(1);
    m = wstar(2);
    E = wstar(3);
    
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    u3 = u2*u;
    p = gm1*(E - rho*u2*0.5);
    h = u*(E+p);
    c = sqrt(gmma*p*orho);
    
    JF_star(3,3) = gmma*u ;
    JF_star(3,2) = gmma*E*orho - gm1*1.5*u2   ;
    JF_star(3,1) = -gmma*u*E*orho + gm1*u3 ;
    JF_star(2,3) = gm1;
    JF_star(2,2) = -gm3*u;
    JF_star(2,1) = gm3*u2*0.5;
    JF_star(1,2) = 1;

    F_star(3,1) = h ;
    F_star(2,1) = m*u + p  ;
    F_star(1,1) = m ;
    
    speedF_star = abs(u)+c;

    dspeeddq_star = 0;
    
    
    
    
    
    speedF_west = speedF_star;    
    truncwest_temp = truncwest_temp - (nuv1*weight)*(vectpsi_west_Trans(:,:,tquad)*F_star);
    if r_param > 0
        
        rho = wwest(1);
        m = wwest(2);
        E = wwest(3);

        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        u3 = u2*u;
        p = gm1*(E - rho*u2*0.5);
        h = u*(E+p);
        c = sqrt(gmma*p*orho);

        JF_west(3,3) = gmma*u ;
        JF_west(3,2) = gmma*E*orho - gm1*1.5*u2   ;
        JF_west(3,1) = -gmma*u*E*orho + gm1*u3 ;
        JF_west(2,3) = gm1;
        JF_west(2,2) = -gm3*u;
        JF_west(2,1) = gm3*u2*0.5;
        JF_west(1,2) = 1;

        F_west(3,1) = h ;
        F_west(2,1) = m*u + p  ;
        F_west(1,1) = m ;

        speedF_west = abs(u)+c;

        dspeeddq_west = 0;        
        wavg = (wwest+wstar)*0.5;
        rho = wavg(1);
        m = wavg(2);
        E = wavg(3);
        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        p = gm1*(E - rho*u2*0.5);        
        c = sqrt(gmma*p*orho);        
        speedF_avg = abs(u)+c;
        
        speedwest = max(speedF_star,speedF_west);
        speedwest = max(speedwest,speedF_avg);
        %speedwest = max([speedF_star speedF_west speedF_avg]);
        Flux = (F_star+F_west-speedwest.*(wstar-wwest))*0.5;
        fluxwest_temp = fluxwest_temp - (nuv1*weight)*(vectpsi_west_Trans(:,:,tquad)*Flux);
    end       
    
    Jac_trunc_west_temp = Jac_trunc_west_temp - vectpsi_west_Trans(:,:,tquad)*JF_star*vectpsi_west(:,:,tquad).*(weight*nuv1);

    if r_param > 0
        %psi_m = vectpsi_west(:,:,tquad);
        dspeeddqr = 0;%compute_F_dwavespeed_dqr(wwest,wstar,quadpoint,d1a);
        F_dqr = (JF_star-speedwest.*I_mat - (wstar-wwest)*dspeeddqr)*0.5;  
        Jac_west_cell_temp = Jac_west_cell_temp - vectpsi_west_Trans(:,:,tquad)*F_dqr*vectpsi_west(:,:,tquad).*weight*nuv1;
        %psi_m = vectpsi_east(:,:,tquad);
        F_dql = 0;%compute_F_dwavespeed_dql(wwest,wstar,quadpoint,d1);
        dspeeddql = F_dql;
        F_dql = (JF_west-speedwest.*I_mat - (wstar-wwest)*dspeeddql)*0.5;  
        Jac_west_other_temp = Jac_west_other_temp - vectpsi_west_Trans(:,:,tquad)*F_dql*vectpsi_east(:,:,tquad).*(weight*nuv1);        
    end
    maxspeedF = max(speedF_star , maxspeedF);
    maxspeedF = max(speedF_west , maxspeedF);
    maxspeedF = max(speedwest1 , maxspeedF);
    maxspeedF = max(speedwest2 , maxspeedF);
end


for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        weight = Dp1quadwgts(k);
        %Derivative terms
                    
        wstar = vectpsi(:,:,tquad,v1quad)*qstar;
        
        rho = wstar(1);
        m = wstar(2);
        E = wstar(3);

        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        u3 = u2*u;
        p = gm1*(E - rho*u2*0.5);
        h = u*(E+p);
        c = sqrt(gmma*p*orho);

        JF_star(3,3) = gmma*u ;
        JF_star(2,3) = gm1;
        JF_star(3,2) = gmma*E*orho - gm1*1.5*u2;
        JF_star(3,1) = -gmma*u*E*orho + gm1*u3 ;
        JF_star(2,2) = -gm3*u;
        JF_star(2,1) = gm3*u2*0.5;
        JF_star(1,2) = 1;

        F_star(3,1) = h ;
        F_star(2,1) = m*u + p  ;
        F_star(1,1) = m ;
        
        speedF_star = abs(u)+c;

        dspeeddq_star = 0;
        
        fterm = (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*F_star;
        residual_cell = residual_cell  - fterm;        
            
        Jac_cell = Jac_cell - (weight*nuv1)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*JF_star*vectpsi(:,:,tquad,v1quad);
        
        maxspeedF = max(speedF_star,maxspeedF);
end






end