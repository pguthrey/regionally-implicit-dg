function [residual,maxspeedF,maxspeedG,maxspeedH] = problem_RKDG_local_residual_eastwest(qstar,qeast,qwest,data,cellcenter)
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
%vectphi_deta_Trans = data.vectphi_deta_Trans;
%vectphi_dzta_Trans = data.vectphi_dzta_Trans;

nuv1 = data.nuv1;

gmma = data.appdata.gamma;
gm1 = gmma-1;
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
            phikdxii = vectphi_dxii_Trans(:,:,v1quad);
            
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
                      
            
            fterm = (nuv1*weight)*phikdxii*Flux;
            system_cell = system_cell - fterm;         
    end
residual = system_cell;
