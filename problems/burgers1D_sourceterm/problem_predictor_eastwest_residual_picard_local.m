function [residual,maxspeedF] ...
          = problem_predictor_eastwest_residual_picard_local(qstar,data,cellcenter)
% Forms the %Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           %Jac_truncdirs
% OUTdata.PUTS   %Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

maxspeedF = 0;

nuv1 = data.nuv1;
Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
Dp1list = data.Dp1list;
Dquadwgts = data.Dquadwgts;
Dp1quadwgts = data.Dp1quadwgts;

vectpsi_east = data.vectpsi_east;
vectpsi_east_Trans = data.vectpsi_east_Trans;
vectpsi_west = data.vectpsi_west;
vectpsi_west_Trans = data.vectpsi_west_Trans;
vectpsi = data.vectpsi;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;

%------------------------------------------------------------------------
%MID Terms

fluxwest_temp = 0;
fluxeast_temp = 0;

for k = 1:Pd
    tquad = Dlist(k,1); 
    %tloc = Dquadlocs(k,1);
    weight = Dquadwgts(k);    

    %East Boundary terms
    wstar = vectpsi_east(:,:,tquad)*qstar;
    
    u = wstar(1);    
    Flux(1,1) = u^2/2 ;    
    speed = abs(u);
    maxspeedF = max(maxspeedF,speed);   
    
    %[speedeast1] = problem_F_Jacobianspectralradius(wstar,1,data.appdata);
    %Flux = problem_F(wstar,1,data.appdata);    
    
    fluxeast_temp = fluxeast_temp + (nuv1*weight)*(vectpsi_east_Trans(:,:,tquad)*Flux);  
    
    %west boundary terms

    wstar = vectpsi_west(:,:,tquad)*qstar;
    
    u = wstar(1);
   
    Flux(1,1) = u^2/2 ;    
    speed = abs(u);
    maxspeedF = max(maxspeedF,speed);   
    
    fluxwest_temp = fluxwest_temp - (nuv1*weight)*(vectpsi_west_Trans(:,:,tquad)*Flux);
end
residual_star = fluxwest_temp + fluxeast_temp;

for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad)*qstar;
        
        u = wstar(1);
        
        Flux(1,1) = u^2/2 ;    
        speed = abs(u);        
        
        residual_star = residual_star - (nuv1*weight)*vectpsi_dxii_Trans(:,:,tquad,v1quad)*Flux;
        maxspeedF = max(speed,maxspeedF);
end

residual =  residual_star ;


end