function [residual,maxspeedF] = problem_RKDG_local_residual_eastwest(qstar,qeast,qwest,data,cellcenter)
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
maxspeedF = 0;

system_cell = zeros(data.theta,1);

vectphi =  data.vectphi;
vectphi_east =  data.vectphi_east;
vectphi_west =  data.vectphi_west;
vectphi_east_Trans =  data.vectphi_east_Trans;
vectphi_west_Trans =  data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;
%vectphi_deta_Trans = data.vectphi_deta_Trans;
%vectphi_dzta_Trans = data.vectphi_dzta_Trans;

nuv1 = data.nuv1;

Pd = data.Pd;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;
Dm1quadwgts = data.Dm1quadwgts;

for k = 1:data.Pdm1
    weight = Dm1quadwgts(k);    

    %East Boundary terms       
    wleft = vectphi_east(:,:)*qstar;
    write = vectphi_west(:,:)*qeast;
    phik = vectphi_east_Trans(:,:);

    u = wleft(1);
    F_left(1,1) = u^2/2;      
    speed = abs(u);

    u = write(1);

    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u));


    wavg = 0.5*(wleft+write);
    u = wavg(1);
    speed = max(speed,abs(u));

    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;    
    speedeast2 = speed;

    %[ Flux ] = problem_F(wstar,quadpoint,data.appdata); 
    %[speedeast1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    %system_cell = system_cell + nuv1*(phik*Flux)*weight;    
    %[ Flux,speedeast2] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
    system_cell = system_cell + (nuv1*weight)*(phik*Flux);  


    %west boundary terms
    wleft = vectphi_east(:,:)*qwest;
    write = vectphi_west(:,:)*qstar;
    phik = vectphi_west_Trans(:,:);

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
    speedwest2 = speed;

    %[Flux] = problem_F(wstar,quadpoint,data.appdata);
    %[speedwest1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
    %system_cell  = system_cell - nuv1*(phik*Flux)*weight;
    %[ Flux,speedwest2 ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
    system_cell = system_cell - (nuv1*weight)*(phik*Flux);

    maxspeedF1 = max(speedeast2,speedwest2);
    maxspeedF = max(maxspeedF,maxspeedF1);
end
for k = 1:Pd
        v1quad = Dlist(k,1); 
        weight = Dquadwgts(k);
        %Derivative terms
        wstar = vectphi(:,:,v1quad)*qstar;

        u = wstar(1);
        Flux(1,1) = u^2/2 ;
        speed = abs(u) ;        
        maxspeedF = max(speed,maxspeedF);

        fterm = (nuv1*weight)*vectphi_dxii_Trans(:,:,v1quad)*Flux;
        system_cell = system_cell - fterm;         
end
residual = system_cell;
