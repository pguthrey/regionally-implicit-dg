function [residual,maxspeedG] = RKDG_local_residual_nortsout(qstar,qnort,qsout,data,cellcenter)
% Forms the Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           Jac_truncdirs
% OUTdata.PUTS   Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    


error('whoops')


maxspeedG = 0;

system_cell = zeros(data.theta,1);

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
        wsout = data.vectphi_nort(:,:,v1quad,v2quad)*qsout;
        wnort = data.vectphi_sout(:,:,v1quad,v2quad)*qnort;

        %East Boundary terms
        v1 = v1center+data.deltav1/2;        
        v2 = v2center+data.deltav2/2*v1loc;         
        v3 = v3center+data.deltav3/2*v2loc;   
        quadpoint = [v1 v2 v3];
        wstar = data.vectphi_nort(:,:,v1quad,v2quad)*qstar;
        phik = data.vectphi_nort_Trans(:,:,v1quad,v2quad);
        %[ Flux ] = problem_F(wstar,quadpoint,data.appdata); 
        %[speednort1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        %system_cell = system_cell + data.nuv2*(phik*Flux)*weight;    
        [ Flux,speednort2] = compute_G_numerical_Flux(wstar,wnort,quadpoint,data.appdata);
        system_cell = system_cell + data.nuv2*(phik*Flux)*weight;  

        speednort1 = 0;
        speedsout1 = 0;

        %sout boundary terms
        v1 = v1center-data.deltav1/2;        
        v2 = v2center+data.deltav2/2*v1loc;         
        v3 = v3center+data.deltav3/2*v2loc;   
        quadpoint = [v1 v2 v3];
        wstar = data.vectphi_sout(:,:,v1quad,v2quad)*qstar;
        phik = data.vectphi_sout_Trans(:,:,v1quad,v2quad);
        %[Flux] = problem_F(wstar,quadpoint,data.appdata);
        %[speedsout1] = problem_F_Jacobianspectralradius(wstar,quadpoint,data.appdata);
        %system_cell  = system_cell - data.nuv2*(phik*Flux)*weight;
        [ Flux,speedsout2 ] = compute_G_numerical_Flux(wsout,wstar,quadpoint,data.appdata);
        system_cell = system_cell - data.nuv2*(phik*Flux)*weight;

        maxspeedG = max([speednort1 speednort2 speedsout1 speedsout2]);
        %wsout = data.vectphi_nort(:,:,v1quad,v2quad)*qsout;
        %wnort = data.vectphi_sout(:,:,v1quad,v2quad)*qnort;
    end

residual = system_cell;
