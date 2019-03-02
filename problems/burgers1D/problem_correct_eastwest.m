function [correct_cell,maxspeedF,data] = problem_correct_eastwest(qstar_1D,qeast_1D,qwest_1D,data,cellcenter,iv1)

correct_cell = 0;
v1quad = 1;

vectpsi_east = data.vectpsi_east;
vectpsi_west = data.vectpsi_west;
vectpsi = data.vectpsi;

vectphi_east_Trans = data.vectphi_east_Trans;
vectphi_west_Trans = data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;

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
    
    u = wleft(1);
    
    F_left(1,1) = u^2/2 ;      
    speed = abs(u) ;
        
    u = write(1);
    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u));        
       
    wavg = 0.5*(wleft+write);
    u = wavg(1);
    speed = max(speed,abs(u));
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;   
    
    maxspeedF = max(maxspeedF,speed); 
    
    correct_cell = correct_cell - (nuv1*weight)*vectphi_east_Trans(:,:,v1quad)*Flux;
    
    %Western Flux
    write = vectpsi_west(:,:,tquad,v1quad)*qstar_1D;
    wleft = vectpsi_east(:,:,tquad,v1quad)*qwest_1D;
        
    u = wleft(1);
    
    F_left(1,1) = u^2/2 ;      
    speed = abs(u) ;
        
    u = write(1);
    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u));        
       
    wavg = 0.5*(wleft+write);
    u = wavg(1);
    speed = max(speed,abs(u) );
            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;       
    maxspeedF = max(maxspeedF,speed);    
   
    correct_cell = correct_cell + (nuv1*weight).*vectphi_west_Trans(:,:,v1quad)*Flux;    
    
end

for k = 1:Pdp1
    tquad = Dp1list(k,1);
    v1quad = Dp1list(k,2);
    weight = Dp1quadwgts(k);
    
    wstar = vectpsi(:,:,tquad,v1quad)*qstar_1D;
    
    u = wstar(1);
    Flux(1,1) = u^2/2 ;
                
    speed = abs(u) ;        
    maxspeedF = max(speed,maxspeedF);        
    
    correct_cell = correct_cell + (nuv1*weight)*vectphi_dxii_Trans(:,:,v1quad)*Flux;
end


end

