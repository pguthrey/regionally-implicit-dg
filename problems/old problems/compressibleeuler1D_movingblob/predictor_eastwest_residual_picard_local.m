function [residual,maxspeedF] ...
          = predictor_eastwest_residual_picard_local(qstar,data,cellcenter)
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

trunceast_temp = 0; 
truncwest_temp = 0; 
fluxeast_temp = 0;
fluxwest_temp = 0;    

maxspeedF = 0;
v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);

speedeast1 = 0;
speedwest1 = 0;
speedeast2 = 0;
speedwest2 = 0;

nuv1 = data.nuv1;
Neqns = data.Neqns;
Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
Dp1list = data.Dp1list;
Dquadwgts = data.Dquadwgts;
Dp1quadlocs = data.Dp1quadlocs;
Dp1quadwgts = data.Dp1quadwgts;
r_param = data.r_param;
space_dims = data.space_dims;

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
    if space_dims >= 2
        v1quad = Dlist(k,2);
        v1loc = Dquadlocs(k,2);
    else
        v1quad = 1;
        v1loc = inf;
    end
    if space_dims >= 3
        v2quad = Dlist(k,3);
        v2loc = Dquadlocs(k,3);
    else
        v2quad = 1;
        v2loc = inf;
    end    
    weight = Dquadwgts(k);    

    %East Boundary terms
    v1 = v1center+deltav1/2;        
    v2 = v2center+deltav2/2*v1loc;         
    v3 = v3center+deltav3/2*v2loc;   
    wstar = vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    psik = vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    [speedeast1] = problem_F_Jacobianspectralradius(wstar,1,appdata);
    Flux = problem_F(wstar,1,appdata);
    fluxeast_temp = fluxeast_temp + nuv1*(psik*Flux)*weight;  
    
    %west boundary terms
    v1 = v1center-deltav1/2;
    v2 = v2center+deltav2/2*v1loc;
    v3 = v3center+deltav3/2*v2loc;

    wstar = vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;
    psik = vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
    [speedwest1] = problem_F_Jacobianspectralradius(wstar,1,appdata);
    Flux = problem_F(wstar,1,appdata);
    fluxwest_temp = fluxwest_temp - nuv1*(psik*Flux)*weight;
    maxspeedF = max([speedeast1 speedwest1 maxspeedF]);
end
residual_star = fluxwest_temp + fluxeast_temp;

v2center = 0;
v3center = 0;

for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        v1loc = Dp1quadlocs(k,2);
        v1center = cellcenter(1);
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            v2loc = Dp1quadlocs(k,3);
            v2center = cellcenter(2);
        else
            v2quad = 1;
            v2loc = inf;
        end
        if space_dims >= 3
            v3quad = Dp1list(k,4);
            v3loc = Dp1quadlocs(k,4);
            v3center = cellcenter(3);
        else
            v3quad = 1;
            v3loc = inf;
        end
        weight = Dp1quadwgts(k);
        %Derivative terms
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        fterm = (nuv1*weight)*psikdxii*problem_F(wstar,1,appdata);
        residual_star = residual_star - fterm;
        [speed2] = problem_F_Jacobianspectralradius(wstar,1,appdata);
        maxspeedF = max([speed2 maxspeedF]);
end

residual =  residual_star ;


end