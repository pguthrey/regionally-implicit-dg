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

vectphi =  data.vectphi;
vectphi_east =  data.vectphi_east;
vectphi_west =  data.vectphi_west;
vectphi_east_Trans =  data.vectphi_east_Trans;
vectphi_west_Trans =  data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;

nuv1 = data.nuv1;

v2quad = 1;
v3quad = 1;
rhs_cell = 0;

Pd = data.Pd;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;

Pdm1 = data.Pdm1;
Dm1list = data.Dm1list;
Dm1quadwgts = data.Dm1quadwgts;

space_dims = data.space_dims;

for k = 1:Pdm1
    if space_dims >= 2
        v1quad = Dm1list(k,1);
        if space_dims >= 3
            v2quad = Dm1list(k,2);
        end
    end
    weight = Dm1quadwgts(k);
    %Flux in the x direction
    %Eastern Flux
    wleft = vectphi_east(:,:,v1quad,v2quad)*qstar;
    write = vectphi_west(:,:,v1quad,v2quad)*qeast;
    phik = vectphi_east_Trans(:,:,v1quad,v2quad);

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

    rhs_cell = rhs_cell + (nuv1*weight)*phik*Flux;

    %Western Flux
    write = vectphi_west(:,:,v1quad,v2quad)*qstar;
    wleft = vectphi_east(:,:,v1quad,v2quad)*qwest;
    phik = vectphi_west_Trans(:,:,v1quad,v2quad);    
    
    u = wleft(1);
    F_left(1,1) = u^2/2 ;      
    speed = abs(u) ;
    u = write(1);     
    F_rite(1,1) = u^2/2 ;   
    speed = max(speed,abs(u) );
    
    Flux = (F_rite+F_left-speed.*(write-wleft))*0.5;    
    speedwest2 = speed;
    
    rhs_cell = rhs_cell - (nuv1*weight).*phik*Flux;
    maxspeedF1 = max(speedeast2,speedwest2);
    maxspeedF = max(maxspeedF,maxspeedF1);    
    
end

for k = 1:Pd
    v1quad = Dlist(k,1);
    if space_dims >= 2
        v2quad = Dlist(k,2);
        if space_dims >= 3
            v3quad = Dlist(k,3);
        end
    end
    weight = Dquadwgts(k);
    wstar = vectphi(:,:,v1quad,v2quad,v3quad)*qstar;
    phikxii = vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
    u = wstar(1);            
    Flux(1,1) = u^2/2 ;
    maxspeedF1 = abs(u);
    maxspeedF = max(maxspeedF,maxspeedF1);        
    rhs_cell = rhs_cell - phikxii*Flux*(nuv1*weight);
end

residual = rhs_cell;
