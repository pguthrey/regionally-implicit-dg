function [residual,maxspeedF] ...
          = problem_predictor_eastwest_residual_picard(qwest,qstar,qeast,data,cellcenter)
% Forms the Jacobian of equations given by the DGmethod for a fixed cell
% written by Pierson Guthrey
% -------------------------------------------------
% INPUTS    DGsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTPUTS   Jacobian
% ------------------------------------------------------------------------    


maxspeedF = 0;


nuv1 = data.nuv1;
%Neqns = data.Neqns;
Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
Dp1list = data.Dp1list;
Dquadwgts = data.Dquadwgts;
%Dp1quadlocs = data.Dp1quadlocs;
Dp1quadwgts = data.Dp1quadwgts;

vectpsi_east = data.vectpsi_east;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west = data.vectpsi_west;
vectpsi_west_Trans = data.vectpsi_west_Trans;
vectpsi = data.vectpsi;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;

gmma = data.appdata.gamma;
gm1 = gmma-1;

%----------------------------------------------------------------
%WEST CELL

fluxeast_temp = 0;
truncwest_temp = 0;

for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    weight = Dquadwgts(k);        
    psik = vectpsi_east_Trans(:,:,tquad);
    %East Boundary terms
    wleft = vectpsi_east(:,:,tquad)*qwest;
    write = vectpsi_west(:,:,tquad)*qstar;
  
    rho = wleft(1);
    m = wleft(2);
    E = wleft(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);   
    h = u*(E+p);    
    F_left(3,1) = h ;
    F_left(2,1) = m*u + p ;
    F_left(1,1) = m ;      
    speed = abs(u) + c;
        
    rho = write(1);
    m = write(2);
    E = write(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    h = u*(E+p);    
    F_rite(3,1) = h ;
    F_rite(2,1) = m*u + p ;
    F_rite(1,1) = m ;   
    speed = max(speed,abs(u) + c);
        
       
    wavg = 0.5*(wleft+write);
    rho = wavg(1);
    m = wavg(2);
    E = wavg(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    speed = max(speed,abs(u) + c);
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;    
    fluxeast_temp = fluxeast_temp + (nuv1*weight)*(psik*Flux);  
    maxspeedF = max(maxspeedF,speed);
    
    %west boundary terms
    wstar = vectpsi_west(:,:,tquad)*qwest;
    psik = vectpsi_west_Trans(:,:,tquad);
    
    rho = wstar(1);
    m = wstar(2);
    E = wstar(3);    
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    h = u*(E+p);
    c = sqrt(gmma*p*orho);    
    Flux(3,1) = h ;
    Flux(2,1) = m*u + p  ;
    Flux(1,1) = m ;    
    speed = abs(u)+c;
    maxspeedF = max(maxspeedF,speed);
    
    truncwest_temp = truncwest_temp - (nuv1*weight)*(psik*Flux);
end

residual_west = fluxeast_temp + truncwest_temp ;

for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2);         
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad)*qwest;
               
        rho = wstar(1);
        m = wstar(2);
        E = wstar(3);

        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        p = gm1*(E - rho*u2*0.5);
        h = u*(E+p);
        c = sqrt(gmma*p*orho);    

        Flux(3,1) = h ;
        Flux(2,1) = m*u + p ;
        Flux(1,1) = m ;
                
        speed = abs(u) + c;        
        maxspeedF = max(speed,maxspeedF);                       
        
        residual_west = residual_west  - (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*Flux;
        %[speed2] = problem_F_Jacobianspectralradius(wstar,1,data.appdata);
end

%------------------------------------------------------------------------
%MID Terms

fluxwest_temp = 0;
fluxeast_temp = 0;

for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    weight = Dquadwgts(k);   
           
    %East Boundary terms    
    wleft = vectpsi_east(:,:,tquad)*qstar;
    write = vectpsi_west(:,:,tquad)*qeast;    
    psik = vectpsi_east_Trans(:,:,tquad);

    rho = wleft(1);
    m = wleft(2);
    E = wleft(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);    
    h = u*(E+p);
    F_left(3,1) = h ;
    F_left(2,1) = m*u + p ;
    F_left(1,1) = m ;      
    speed = abs(u) + c;
        
    rho = write(1);
    m = write(2);
    E = write(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    h = u*(E+p);
    F_rite(3,1) = h ;
    F_rite(2,1) = m*u + p ;
    F_rite(1,1) = m ;   
    speed = max(speed,abs(u) + c);
        
       
    wavg = 0.5*(wleft+write);
    rho = wavg(1);
    m = wavg(2);
    E = wavg(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    speed = max(speed,abs(u) + c);
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;    
    fluxeast_temp = fluxeast_temp + (nuv1*weight)*(psik*Flux);  
    maxspeedF = max(maxspeedF,speed);
    
    %west boundary terms
    wleft = vectpsi_east(:,:,tquad)*qwest;
    write = vectpsi_west(:,:,tquad)*qstar;
    psik = vectpsi_west_Trans(:,:,tquad);
        
    rho = wleft(1);
    m = wleft(2);
    E = wleft(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);    
    h = u*(E+p);
    F_left(3,1) = h ;
    F_left(2,1) = m*u + p ;
    F_left(1,1) = m ;      
    speed = abs(u) + c;
        
    rho = write(1);
    m = write(2);
    E = write(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    F_rite(3,1) = h ;
    F_rite(2,1) = m*u + p ;
    F_rite(1,1) = m ;   
    speed = max(speed,abs(u) + c);
        
       
    wavg = 0.5*(wleft+write);
    rho = wavg(1);
    m = wavg(2);
    E = wavg(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    speed = max(speed,abs(u) + c);
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;                 
    fluxwest_temp = fluxwest_temp - (nuv1*weight)*(psik*Flux);
    
    maxspeedF = max(maxspeedF,speed);    
end
residual_star = fluxwest_temp + fluxeast_temp;

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
        p = gm1*(E - rho*u2*0.5);
        c = sqrt(gmma*p*orho);    
        h = u*(E+p);
        Flux(3,1) = h ;
        Flux(2,1) = m*u + p ;
        Flux(1,1) = m ;
                
        speed = abs(u) + c;        
        maxspeedF = max(speed,maxspeedF);
        
        residual_star = residual_star - (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*Flux;
        %[speed2] = problem_F_Jacobianspectralradius(wstar,1,data.appdata);
        
end

%------------------------------------------------------------------------
%EAST Terms

trunceast_temp = 0;
fluxwest_temp = 0;

for k = 1:Pd
    tquad = Dlist(k,1); 
    weight = Dquadwgts(k);    
    
    %East Boundary terms
    wstar = vectpsi_east(:,:,tquad)*qeast;
    psik = vectpsi_east_Trans(:,:,tquad);

    rho = wstar(1);
    m = wstar(2);
    E = wstar(3);    
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);    
    h = u*(E+p);
    Flux(3,1) = h ;
    Flux(2,1) = m*u + p  ;
    Flux(1,1) = m ;    
    speed = abs(u)+c;
    maxspeedF = max(maxspeedF,speed);
    
    trunceast_temp = trunceast_temp + (nuv1*weight)*(psik*Flux);
    
    %west boundary terms
    wleft = vectpsi_east(:,:,tquad)*qstar;
    write = vectpsi_west(:,:,tquad)*qeast;   
    psik = vectpsi_west_Trans(:,:,tquad);
    
    rho = wleft(1);
    m = wleft(2);
    E = wleft(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);    
    h = u*(E+p);
    F_left(3,1) = h ;
    F_left(2,1) = m*u + p ;
    F_left(1,1) = m ;      
    speed = abs(u) + c;
        
    rho = write(1);
    m = write(2);
    E = write(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    h = u*(E+p);
    F_rite(3,1) = h ;
    F_rite(2,1) = m*u + p ;
    F_rite(1,1) = m ;   
    speed = max(speed,abs(u) + c);
        
       
    wavg = 0.5*(wleft+write);
    rho = wavg(1);
    m = wavg(2);
    E = wavg(3);
    orho = 1/rho;
    u = m*orho;
    u2 = u*u;
    p = gm1*(E - rho*u2*0.5);
    c = sqrt(gmma*p*orho);
    speed = max(speed,abs(u) + c);
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;         
    
    fluxwest_temp = fluxwest_temp - (nuv1*weight)*(psik*Flux);
    
    maxspeedF = max(maxspeedF,speed);
end
residual_east = trunceast_temp + fluxwest_temp;

%%{
for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad)*qeast;
        
        rho = wstar(1);
        m = wstar(2);
        E = wstar(3);

        orho = 1/rho;
        u = m*orho;
        u2 = u*u;
        p = gm1*(E - rho*u2*0.5);
        c = sqrt(gmma*p*orho);    
        h = u*(E+p);

        Flux(3,1) = h ;
        Flux(2,1) = m*u + p ;
        Flux(1,1) = m ;
                
        speed = abs(u) + c;        
        maxspeedF = max(speed,maxspeedF);       
                
        residual_east = residual_east  - (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*Flux;
end
%}

residual = [residual_west ; residual_star ; residual_east ];



end