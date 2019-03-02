function [residual,maxspeedG] = problem_RKDG_local_residual_nortsout(qstar,qnort,qsout,data,cellcenter)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INPUTS    DGsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTPUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

%quickfix
maxspeedH = 0;
maxspeedF = 0;
maxspeedG = 0;

system_cell = zeros(data.theta,1);

vectphi =  data.vectphi;
vectphi_east =  data.vectphi_east;
vectphi_west =  data.vectphi_west;
vectphi_east_Trans =  data.vectphi_east_Trans;
vectphi_west_Trans =  data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;

vectphi_nort =  data.vectphi_nort;
vectphi_sout =  data.vectphi_sout;
vectphi_nort_Trans =  data.vectphi_nort_Trans;
vectphi_sout_Trans =  data.vectphi_sout_Trans;
vectphi_deta_Trans = data.vectphi_deta_Trans;

nuv1 = data.nuv1;
nuv2 = data.nuv2;

Pd = data.Pd;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;
Dm1quadwgts = data.Dm1quadwgts;

correct_cell = 0;

% North - south contribution 
%----------------------------
for k = 1:data.Pd
    weight = data.Dquadwgts(k);
    v1quad = data.Dlist(k,1);    
    v2quad = data.Dlist(k,2);    
    %Flux in the y direction
    %Northern Flux
    wleft = vectphi_nort(:,:,v1quad,v2quad)*qstar;
    write = vectphi_sout(:,:,v1quad,v2quad)*qnort;
    phik = vectphi_nort_Trans(:,:,v1quad,v2quad);

    u = wleft(1);
    F_left(1,1) = u^2/2 ;      
    speed = abs(u) ;
    u = write(1);     
    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u) );
    wavg = 0.5*(wleft+write);
    u = wavg(1);        
    speed = max(speed,abs(u) );            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;    
        
    correct_cell = correct_cell + (nuv2*weight)*phik*Flux;
            
    speednort = speed;
    
    %Southern Flux
    write = data.vectphi_sout(:,:,v1quad,v2quad)*qstar;
    wleft = data.vectphi_nort(:,:,v1quad,v2quad)*qsout;
    phik = data.vectphi_sout_Trans(:,:,v1quad,v2quad);
    
    u = wleft(1);
    F_left(1,1) = u^2/2 ;      
    speed = abs(u) ;
    u = write(1);     
    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u) );
    wavg = 0.5*(wleft+write);
    u = wavg(1);        
    speed = max(speed,abs(u) );            
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;  
        
    correct_cell = correct_cell - (nuv2*weight)*phik*Flux;            
    
    speedsout = speed;
    
    maxspeedG = max([speednort speedsout]);

end
% derivative contribution
for k = 1:Pd
        v1quad = Dlist(k,1); 
        v2quad = Dlist(k,2); 
        weight = Dquadwgts(k);
        %Derivative terms
        wstar = vectphi(:,:,v1quad,v2quad)*qstar;
        %phikdxii = ;

        u = wstar(1);            
        Flux(1,1) = u^2/2 ;

        speed = abs(u) ;        
        maxspeedF = max(speed,maxspeedF);                   

        fterm = (nuv1*weight)*vectphi_dxii_Trans(:,:,v1quad)*Flux;
        system_cell = system_cell - fterm;         
end

residual = correct_cell;