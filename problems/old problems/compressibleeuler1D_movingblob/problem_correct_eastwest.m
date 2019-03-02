function [correct_cell,maxspeedF,data] = problem_correct_eastwest(qstar_1D,qeast_1D,qwest_1D,data,cellcenter,iv1)

correct_cell = 0;
v1quad = 1;

vectpsi_east = data.vectpsi_east;
vectpsi_west = data.vectpsi_west;
vectpsi = data.vectpsi;

vectphi_east_Trans = data.vectphi_east_Trans;
vectphi_west_Trans = data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;


gmma = data.appdata.gamma;
gm1 = gmma-1;

nuv1 = data.nuv1;
%Neqns = data.Neqns;
Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
Dp1list = data.Dp1list;
Dquadwgts = data.Dquadwgts;
%Dp1quadlocs = data.Dp1quadlocs;
Dp1quadwgts = data.Dp1quadwgts;

maxspeedF = 0;

for k = 1:Pd
    weight = Dquadwgts(k);
    tquad = Dlist(k,1);
    %Flux in the x direction
    %Eastern Flux
    wleft = vectpsi_east(:,:,tquad,v1quad)*qstar_1D;
    write = vectpsi_west(:,:,tquad,v1quad)*qeast_1D;
    
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
    
    maxspeedF = max(maxspeedF,speed); 
    
    correct_cell = correct_cell - (nuv1*weight)*vectphi_east_Trans(:,:,v1quad)*Flux;
    
    %Western Flux
    write = vectpsi_west(:,:,tquad,v1quad)*qstar_1D;
    wleft = vectpsi_east(:,:,tquad,v1quad)*qwest_1D;
        
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
    maxspeedF = max(maxspeedF,speed);    
   
    correct_cell = correct_cell + (nuv1*weight).*vectphi_west_Trans(:,:,v1quad)*Flux;    
    
end

for k = 1:Pdp1
    tquad = Dp1list(k,1);
    v1quad = Dp1list(k,2);
    weight = Dp1quadwgts(k);
    
    wstar = vectpsi(:,:,tquad,v1quad)*qstar_1D;
    
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
    
    correct_cell = correct_cell + (nuv1*weight)*vectphi_dxii_Trans(:,:,v1quad)*Flux;
end
end

