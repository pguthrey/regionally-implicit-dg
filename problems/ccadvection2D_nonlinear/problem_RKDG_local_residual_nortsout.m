function [residual,maxspeedG] = problem_RKDG_local_residual_nortsout(qstar,qnort,qsout,correct_cell,data,cellcenter)
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

maxspeedG = 0;

vectphi =  data.vectphi;
vectphi_nort =  data.vectphi_nort;
vectphi_sout =  data.vectphi_sout;
vectphi_nort_Trans =  data.vectphi_nort_Trans;
vectphi_sout_Trans =  data.vectphi_sout_Trans;
vectphi_deta_Trans = data.vectphi_deta_Trans;

nuv2 = data.nuv2;

v2quad = 1;
v3quad = 1;

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
    %Glux in the x direction
    %Eastern Glux
    wleft = vectphi_nort(:,:,v1quad,v2quad)*qstar;
    write = vectphi_sout(:,:,v1quad,v2quad)*qnort;
    phik = vectphi_nort_Trans(:,:,v1quad,v2quad);

    u = wleft(1);
    G_left(1,1) = data.appdata.nug*u;      
    speed = abs(data.appdata.nug);

    u = write(1);

    G_rite(1,1) = data.appdata.nug*u;   
    speed = max(speed,abs(data.appdata.nug));

    wavg = 0.5*(wleft+write);
    u = wavg(1);
    speed = max(speed,abs(data.appdata.nug));

    Glux = (G_rite+G_left-speed.*(write-wleft))*0.5;    
    speednort2 = speed;

    correct_cell = correct_cell + (nuv2*weight)*phik*Glux;

    %Western Glux
    write = vectphi_sout(:,:,v1quad,v2quad)*qstar;
    wleft = vectphi_nort(:,:,v1quad,v2quad)*qsout;
    phik = vectphi_sout_Trans(:,:,v1quad,v2quad);    
    
    u = wleft(1);
    G_left(1,1) = data.appdata.nug*u;      
    speed = abs(data.appdata.nug) ;
    u = write(1);     
    G_rite(1,1) = data.appdata.nug*u ;   
    speed = max(speed,abs(data.appdata.nug));
    
    Glux = (G_rite+G_left-speed.*(write-wleft))*0.5;    
    speedsout2 = speed;
    
    correct_cell = correct_cell - (nuv2*weight).*phik*Glux;
    maxspeedG1 = max(speednort2,speedsout2);
    maxspeedG = max(maxspeedG,maxspeedG1);    
    
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
    phiketa = vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
    u = wstar(1);            
    Glux(1,1) = data.appdata.nug*u;
    correct_cell = correct_cell - phiketa*Glux*(nuv2*weight);
end

residual = correct_cell;
