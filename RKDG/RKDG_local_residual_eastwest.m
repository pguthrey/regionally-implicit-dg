function [residual,maxspeedF,maxspeedG,maxspeedH] = RKDG_local_residual_eastwest(qstar,qeast,qwest,data,cellcenter)
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


error('whoops')
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
%vectphi_deta_Trans = data.vectphi_deta_Trans;
%vectphi_dzta_Trans = data.vectphi_dzta_Trans;

nuv1 = data.nuv1;

v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);
   
    for k = 1:data.Pdm1
        if data.space_dims >= 2
            v1quad = data.Dm1list(k,1);
            v1loc = data.Dm1quadlocs(k,1);
        else
            v1quad = 1;
            v1loc = inf;
        end
        if data.space_dims >= 3
            v2quad = data.Dm1list(k,2);
            v2loc = data.Dm1quadlocs(k,2);
        else
            v2quad = 1;
            v2loc = inf;
        end    
        weight = data.Dm1quadwgts(k);    
        wwest = vectphi_east(:,:,v1quad,v2quad)*qwest;
        weast = vectphi_west(:,:,v1quad,v2quad)*qeast;

        %East Boundary terms
        v1 = v1center+data.deltav1/2;        
        v2 = v2center+data.deltav2/2*v1loc;         
        v3 = v3center+data.deltav3/2*v2loc;   
        quadpoint = [v1 v2 v3];
        wstar = vectphi_east(:,:,v1quad,v2quad)*qstar;
        phik = vectphi_east_Trans(:,:,v1quad,v2quad);
        %[ Flux ] = problem_F(wstar,quadpoint,data.appdata); 
        %[speedeast1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        %system_cell = system_cell + nuv1*(phik*Flux)*weight;    
        [ Flux,speedeast2] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
        system_cell = system_cell + (nuv1*weight)*(phik*Flux);  

        speedeast1 = 0;
        speedwest1 = 0;

        %west boundary terms
        v1 = v1center-data.deltav1/2;        
        v2 = v2center+data.deltav2/2*v1loc;         
        v3 = v3center+data.deltav3/2*v2loc;   
        quadpoint = [v1 v2 v3];
        wstar = vectphi_west(:,:,v1quad,v2quad)*qstar;
        phik = vectphi_west_Trans(:,:,v1quad,v2quad);
        %[Flux] = problem_F(wstar,quadpoint,data.appdata);
        %[speedwest1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        %system_cell  = system_cell - nuv1*(phik*Flux)*weight;
        [ Flux,speedwest2 ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
        system_cell = system_cell - (nuv1*weight)*(phik*Flux);

        maxspeedF = max([speedeast1 speedeast2 speedwest1 speedwest2]);
    end

    v2center = 0;
    v3center = 0;

    for k = 1:data.Pd
            v1quad = data.Dlist(k,1); 
            v1loc = data.Dquadlocs(k,1);
            v1center = cellcenter(1);
            if data.space_dims >= 2
                v2quad = data.Dlist(k,2);
                v2loc = data.Dquadlocs(k,2);
                v2center = cellcenter(2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.Dlist(k,3);
                v3loc = data.Dquadlocs(k,3);
                v3center = cellcenter(3);
            else
                v3quad = 1;
                v3loc = inf;
            end
            weight = data.Dquadwgts(k);
            %Derivative terms
            v1 = v1center+data.deltav1/2*v1loc;
            v2 = v2center+data.deltav2/2*v2loc;
            v3 = v3center+data.deltav3/2*v3loc;
            quadpoint = [v1 v2 v3];
            wstar = vectphi(:,:,v1quad,v2quad,v3quad)*qstar;
            phikdxii = vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
            fterm = (nuv1*weight)*phikdxii*problem_F(wstar,quadpoint,data.appdata);
            system_cell = system_cell - fterm;         

    end
residual = system_cell;
