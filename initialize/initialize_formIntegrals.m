function [data] = initialize_formIntegrals(data,timeleft)    
% Pre-evaluates the integrals used in the tau integration-by-parts 
% written by Pierson Guthrey

if ~strcmp(data.linearcase,'constantcoefficient') 
    data.I_futr = zeros(data.thetaT,data.thetaT) ;
    data.I_past = zeros(data.thetaT,data.thetaT);
    data.Psi_dtau = zeros(data.thetaT,data.thetaT);

    data.I_futr = zeros(data.thetaT,data.thetaT);                
    data.I_past = zeros(data.thetaT,data.thetaT); 
    data.Psi_dtau = zeros(data.thetaT,data.thetaT);                            

    temp1 = 0;
    temp2 = 0;
    temp3 = 0;

    for k = 1:data.Pd
        v1quad = data.Dlist(k,1); 
        if data.space_dims >= 2
            v2quad = data.Dlist(k,2);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
        else
            v3quad = 1;
        end
        weight = data.Dquadwgts(k);
        psikplus = data.vectpsi_futr_Trans(:,:,v1quad,v2quad,v3quad);
        psikmnus = data.vectpsi_past_Trans(:,:,v1quad,v2quad,v3quad);
        psil = data.vectpsi_futr(:,:,v1quad,v2quad,v3quad);
        temp1 = temp1 + psikplus*psil*weight;
        temp2 = temp2 + psikmnus*psil*weight;
    end

    for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        v1quad = data.Dp1list(k,2); 
        if data.space_dims >= 2
            v2quad = data.Dp1list(k,3);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dp1list(k,4);
        else
            v3quad = 1;
        end
        weight = data.Dp1quadwgts(k);
        psil = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        psiktau = data.vectpsi_dtau_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        temp3 = temp3 + psiktau*psil*weight;
    end        
    data.I_futr = temp1;                
    data.I_past = temp2;
    data.Psi_dtau = temp3;
    data.tauflux = data.I_futr-data.Psi_dtau;                


    projection_temp = 0;
    PHI_norm_temp = 0;
    PSI_norm_temp = 0;

    for v1quad = 1:data.n1quad
    for v2quad = 1:data.n2quad
    for v3quad = 1:data.n3quad
        for tquad = 1:data.n1quad
            weight = data.wgts_spacetime(tquad,v1quad,v2quad,v3quad);
            psitrans = data.vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad);
            psi = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
            phi = data.vectphi(:,:,v1quad,v2quad,v3quad);

            projection_temp = projection_temp + psitrans*phi*weight;
            PSI_norm_temp = PSI_norm_temp + psitrans*psi*weight;
        end
    end
    end
    end

    for v1quad = 1:data.n1quad
    for v2quad = 1:data.n2quad
    for v3quad = 1:data.n3quad
            weight = data.wgts_space(v1quad,v2quad,v3quad);
            phitrans = data.vectphi_Trans(:,:,v1quad,v2quad,v3quad);
            phi = data.vectphi(:,:,v1quad,v2quad,v3quad);
            PHI_norm_temp = PHI_norm_temp + phitrans*phi*weight;
    end
    end
    end

    data.PHI_norm = PHI_norm_temp;
    data.PSI_norm = PSI_norm_temp;
    data.PHI_norm_inv = inv(PHI_norm_temp);
    data.PSI_norm_inv = inv(PSI_norm_temp);
    data.PHI2PSI = data.PSI_norm_inv*projection_temp;


    data.extrapolate_psi_x = 0;
    data.extrapolate_psi_y = 0;
    for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        v1quad = data.Dp1list(k,2); 
        p = max(data.Dp1list(:,2));
        if data.space_dims >= 2
            v2quad = data.Dp1list(k,3);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dp1list(k,4);
        else
            v3quad = 1;
        end    
        weight = data.Dp1quadwgts(k);
        psiT = data.vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        psi = data.vectpsi(:,:,tquad,p+1-v1quad,v2quad,v3quad);
        data.extrapolate_psi_x = data.extrapolate_psi_x + weight*psiT*psi;
        if data.space_dims >= 2
            psi = data.vectpsi(:,:,tquad,v1quad,data.P+1-v2quad,v3quad);
            data.extrapolate_psi_y = data.extrapolate_psi_y + weight*psiT*psi;
        end
    end
    data.extrapolate_psi_x = data.PSI_norm_inv*data.extrapolate_psi_x;
    if data.space_dims >= 2
        data.extrapolate_psi_y = data.PSI_norm_inv*data.extrapolate_psi_y;
    end

    %{
    data.extrapolate_phi_x = 0;
    data.extrapolate_phi_y = 0;
    for k = 1:data.Pd
        v1quad = data.Dlist(k,1); 
        v1loc = data.Dquadlocs(k,1);
        if data.space_dims >= 2
            v2quad = data.Dlist(k,2);
            v2loc = data.Dquadlocs(k,2);
        else
            v2quad = 1;
            v2loc = inf;
        end
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
            v3loc = data.Dquadlocs(k,3);
        else
            v3quad = 1;
            v3loc = inf;
        end    
        quadpoint = [v1loc,v2loc,v3loc];
        weight = data.Dquadwgts(k);
        phiT = data.vectphi_Trans(:,:,v1quad,v2quad,v3quad);
        phi = data.vectphi(:,:,data.P+1-v1quad,v2quad,v3quad);
        data.extrapolate_phi_x = data.extrapolate_phi_x + weight*phiT*phi;
        if data.space_dims >= 2
            phi = data.vectphi(:,:,v1quad,data.P+1-v2quad,v3quad);
            data.extrapolate_phi_y = data.extrapolate_phi_y + weight*phiT*phi;
        end
    end
    data.extrapolate_phi_x = data.PHI_norm_inv*data.extrapolate_phi_x;
    if data.space_dims >= 2
        data.extrapolate_phi_y = data.PHI_norm_inv*data.extrapolate_phi_y;
    end
    %}


    if 0%data.linearcase
        r = data.r_param;
        ry = 0;
        rz = 0;
        if data.space_dims >= 2
            ry = data.r_param;
            if data.space_dims >= 3
                rz = data.r_param;
            end
        end

        index = @(gamma) (gamma-1)*data.thetaT+(1:data.thetaT);
        regions = 2*data.r_param+1;
        matsize = regions^data.space_dims*data.thetaT;
        data.update = NaN(data.thetaT,matsize,data.Nv1,data.Nv2,data.Nv3);
        globalspeedF = 0;               
        globalspeedG = 0;
        globalspeedH = 0;


        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2
        for iv3 = 1:data.Nv3
            for ir3 = -rz:rz            
            for ir2 = -ry:ry
            for ir1 = -r:r
                for k = 1:data.Pd
                    if data.space_dims >= 2
                        v1loc = data.Dquadlocs(k,2);
                    else
                        v1loc = inf;
                    end
                    if data.space_dims >= 3
                        v2loc = data.Dquadlocs(k,3);
                    else
                        v2loc = inf;
                    end    

                    %East Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                    quadpoint = [v1 v2 v3];

                    localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    globalspeedF = max([localspeedF globalspeedF]);

                    %West Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1-data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                    quadpoint = [v1 v2 v3];

                    localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    globalspeedF = max([localspeedF globalspeedF]);

                    if data.space_dims >= 2
                        %nort Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                        quadpoint = [v1 v2 v3];

                        localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        globalspeedG = max([localspeedG globalspeedG]);

                        %sout Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2-data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                        quadpoint = [v1 v2 v3];

                        localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        globalspeedG = max([localspeedG globalspeedG]);
                    end
                end


                for k = 1:data.Pdp1
                    %tloc = data.Dp1quadlocs(k,1);
                    v1loc = data.Dp1quadlocs(k,2);
                    if data.space_dims >= 2
                        v2loc = data.Dp1quadlocs(k,3);
                    else
                        v2loc = inf;
                    end
                    if data.space_dims >= 3
                        v3loc = data.Dp1quadlocs(k,4);
                    else
                        v3loc = inf;
                    end
                    %Derivative terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v2loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v3loc;   
                    quadpoint = [v1 v2 v3];

                    localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    globalspeedF = max([localspeedF globalspeedF]);
                    if data.space_dims >= 2
                        localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        globalspeedG = max([localspeedG globalspeedG]);
                    end     
                end             
            end
            end
            end
        end
        end
        end
        data.Fspeedmax = globalspeedF;
        data.Gspeedmax = globalspeedG;
        data.Hspeedmax = globalspeedH;

        deltat1 = data.deltav1/data.Fspeedmax*data.cfl;
        deltat2 = data.deltav2/data.Gspeedmax*data.cfl;
        deltat3 = data.deltav3/data.Hspeedmax*data.cfl;
        deltat = min([deltat1 deltat2 deltat3 timeleft]);
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;
        data.nuv3 = deltat/data.deltav3;

        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2
        for iv3 = 1:data.Nv3
            cellnum = 1;
            LHS = zeros(matsize);
            RHS = LHS;
            for ir3 = -rz:rz            
            for ir2 = -ry:ry
            for ir1 = -r:r
                row = index(cellnum);
                col = index(cellnum);

                fluxeast_temp = 0;                
                trunceast_temp = 0;
                othereast_temp = 0;
                fluxwest_temp = 0;                
                truncwest_temp = 0;  
                otherwest_temp = 0;
                jaccell_temp = 0;

                fluxnort_temp = 0;                
                truncnort_temp = 0;
                othernort_temp = 0;
                fluxsout_temp = 0;                
                truncsout_temp = 0;  
                othersout_temp = 0;

                %{
                fluxuppr_temp = 0;                
                truncuppr_temp = 0;
                otheruppr_temp = 0;
                fluxdown_temp = 0;                
                truncdown_temp = 0;  
                otherdown_temp = 0;
                %}


                for k = 1:data.Pd
                    tquad = data.Dlist(k,1); 
                    %tloc = data.Dquadlocs(k,1);
                    if data.space_dims >= 2
                        v1quad = data.Dlist(k,2);
                        v1loc = data.Dquadlocs(k,2);
                    else
                        v1quad = 1;
                        v1loc = inf;
                    end
                    if data.space_dims >= 3
                        v2quad = data.Dlist(k,3);
                        v2loc = data.Dquadlocs(k,3);
                    else
                        v2quad = 1;
                        v2loc = inf;
                    end    
                    weight = data.Dquadwgts(k);

                    %East Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                    quadpoint = [v1 v2 v3];

                    psik = data.vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
                    psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    Flux = predictor_F_dql(NaN,NaN,quadpoint,data.appdata); 
                    fluxeast_temp = fluxeast_temp + data.nuv1*psik*Flux*psil*weight;                
                    Flux = f;
                    trunceast_temp = trunceast_temp + data.nuv1*psik*Flux*psil*weight;
                    psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
                    Flux = predictor_F_dqr(NaN,NaN,quadpoint,data.appdata);
                    othereast_temp = othereast_temp + data.nuv1*psik*Flux*psil*weight;                
                    globalspeedF = max([localspeedF globalspeedF]);

                    %West Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1-data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                    quadpoint = [v1 v2 v3];

                    psik = data.vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
                    psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    Flux = predictor_F_dqr(NaN,NaN,quadpoint,data.appdata);
                    fluxwest_temp = fluxwest_temp + data.nuv1*psik*Flux*psil*weight;                
                    Flux = f;
                    truncwest_temp = truncwest_temp + data.nuv1*psik*Flux*psil*weight;
                    psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
                    Flux = predictor_F_dql(NaN,NaN,quadpoint,data.appdata);
                    otherwest_temp = otherwest_temp + data.nuv1*psik*Flux*psil*weight;
                    globalspeedF = max([localspeedF globalspeedF]);

                    if data.space_dims >= 2
                        %nort Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                        quadpoint = [v1 v2 v3];

                        psik = data.vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
                        psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        Flux = predictor_G_dql(NaN,NaN,quadpoint,data.appdata);
                        fluxnort_temp = fluxnort_temp + data.nuv2*psik*Flux*psil*weight;                
                        Flux = g;
                        truncnort_temp = truncnort_temp + data.nuv2*psik*Flux*psil*weight;
                        psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                        Flux = predictor_G_dqr(NaN,NaN,quadpoint,data.appdata);
                        othernort_temp = othernort_temp + data.nuv2*psik*Flux*psil*weight;                
                        globalspeedG = max([localspeedG globalspeedG]);

                        %sout Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2-data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                        quadpoint = [v1 v2 v3];

                        psik = data.vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
                        psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        Flux = predictor_G_dqr(NaN,NaN,quadpoint,data.appdata);
                        fluxsout_temp = fluxsout_temp + data.nuv2*psik*Flux*psil*weight;                
                        Flux = g;
                        truncsout_temp = truncsout_temp + data.nuv2*psik*Flux*psil*weight;
                        psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                        Flux = predictor_G_dql(NaN,NaN,quadpoint,data.appdata);
                        othersout_temp = othersout_temp + data.nuv2*psik*Flux*psil*weight;
                        globalspeedG = max([localspeedG globalspeedG]);
                    end
                end


                for k = 1:data.Pdp1
                    tquad = data.Dp1list(k,1); 
                    %tloc = data.Dp1quadlocs(k,1);
                    v1quad = data.Dp1list(k,2); 
                    v1loc = data.Dp1quadlocs(k,2);
                    if data.space_dims >= 2
                        v2quad = data.Dp1list(k,3);
                        v2loc = data.Dp1quadlocs(k,3);
                    else
                        v2quad = 1;
                        v2loc = inf;
                    end
                    if data.space_dims >= 3
                        v3quad = data.Dp1list(k,4);
                        v3loc = data.Dp1quadlocs(k,4);
                    else
                        v3quad = 1;
                        v3loc = inf;
                    end
                    weight = data.Dp1quadwgts(k);
                    %Derivative terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v2loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v3loc;   
                    quadpoint = [v1 v2 v3];

                    phil = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
                    psikdxii = data.vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    globalspeedF = max([localspeedF globalspeedF]);
                    fterm = weight*data.nuv1*psikdxii*f*phil;
                    jaccell_temp = jaccell_temp + fterm;
                    if data.space_dims >= 2
                        psikdeta = data.vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        globalspeedG = max([localspeedG globalspeedG]);
                        gterm = weight*data.nuv2*psikdeta*g*phil;
                        jaccell_temp = jaccell_temp + gterm;
                    end               
                end


                LHS(row,col) = data.I_futr - data.Psi_dtau -  jaccell_temp;
                RHS(row,col) = data.I_past;                
                %West LHS flux
                if ir1 > -r
                    LHS(row,col) = LHS(row,col) - fluxwest_temp;
                    col = index(cellnum-1);
                    LHS(row,col) = LHS(row,col) - otherwest_temp;
                else
                    LHS(row,col) = LHS(row,col) - truncwest_temp;
                end
                col = index(cellnum);
                %East LHS flux
                if ir1 < r
                    LHS(row,col) = LHS(row,col) + fluxeast_temp;
                    col = index(cellnum+1);
                    LHS(row,col) = LHS(row,col) + othereast_temp;
                else
                    LHS(row,col) = LHS(row,col) + trunceast_temp;
                end            
                col = index(cellnum);

                if data.space_dims >= 2
                    %sout LHS flux
                    if ir2 > -ry
                        LHS(row,col) = LHS(row,col) - fluxsout_temp;
                        col = index(cellnum-regions);
                        LHS(row,col) = LHS(row,col) - othersout_temp;
                    else
                        LHS(row,col) = LHS(row,col) - truncsout_temp;
                    end
                    col = index(cellnum);
                    %nort LHS flux
                    if ir2 < ry
                        LHS(row,col) = LHS(row,col) + fluxnort_temp;
                        col = index(cellnum+regions);
                        LHS(row,col) = LHS(row,col) + othernort_temp;
                    else
                        LHS(row,col) = LHS(row,col) + truncnort_temp;
                    end            
                    if data.space_dims >= 3
                        col = index(cellnum);
                        %down LHS flux
                        if ir3 > -rz
                            LHS(row,col) = LHS(row,col) + this.Jac_down_cell(:,:,celli,cellj,cellk);
                            col = index(cellnum-regions^2);
                            LHS(row,col) = LHS(row,col) + this.Jac_down_other(:,:,celli,cellj,cellk);
                        else
                            LHS(row,col) = LHS(row,col) + this.Jac_trunc_down(:,:,celli,cellj,cellk);
                        end
                        col = index(cellnum);
                        %uppr LHS flux
                        if ir3 < rz
                            LHS(row,col) = LHS(row,col) + this.Jac_uppr_cell(:,:,celli,cellj,cellk);
                            col = index(cellnum+regions^2);
                            LHS(row,col) = LHS(row,col) + this.Jac_uppr_other(:,:,celli,cellj,cellk);
                        else
                            LHS(row,col) = LHS(row,col) + this.Jac_trunc_uppr(:,:,celli,cellj,cellk);
                        end  
                    end
                end
                %iterate cell number
                cellnum = cellnum + 1;                
            end
            end
            end
            coeffs = (1:data.thetaT)+(data.main_cell-1)*data.thetaT;

            if iv1 == 24
            %    keyboard
            end

            solution = LHS\RHS;
            %keyboard
            data.predictor_update(:,:,iv1,iv2,iv3) = solution(coeffs,:);        

        end
        end
        end

        data.corrector_update =  zeros(data.theta,data.thetaT*(1+2*data.space_dims),data.Nv1,data.Nv2,data.Nv3);
        globalspeedF = 0;
        globalspeedG = 0;

        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2
        for iv3 = 1:data.Nv3

                fluxeast_temp = 0;       
                othereast_temp = 0;
                fluxwest_temp = 0;                
                otherwest_temp = 0;
                jaccell_temp = 0;                      

                fluxnort_temp = 0;                
                othernort_temp = 0;
                fluxsout_temp = 0;                
                othersout_temp = 0;

                %{
                fluxuppr_temp = 0;                
                truncuppr_temp = 0;
                otheruppr_temp = 0;
                fluxdown_temp = 0;                
                truncdown_temp = 0;  
                otherdown_temp = 0;
                %}


                for k = 1:data.Pd
                    tquad = data.Dlist(k,1); 
                    %tloc = data.Dquadlocs(k,1);
                    if data.space_dims >= 2
                        v1quad = data.Dlist(k,2);
                        v1loc = data.Dquadlocs(k,2);
                    else
                        v1quad = 1;
                        v1loc = inf;
                    end
                    if data.space_dims >= 3
                        v2quad = data.Dlist(k,3);
                        v2loc = data.Dquadlocs(k,3);
                    else
                        v2quad = 1;
                        v2loc = inf;
                    end    
                    weight = data.Dquadwgts(k);

                    %East Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                    quadpoint = [v1 v2 v3];

                    psik = data.vectphi_east_Trans(:,:,v1quad,v2quad);
                    psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    Flux = predictor_F_dql(NaN,NaN,quadpoint,data.appdata);
                    fluxeast_temp = fluxeast_temp + data.nuv1*psik*Flux*psil*weight;                
                    psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
                    Flux = predictor_F_dqr(NaN,NaN,quadpoint,data.appdata);
                    othereast_temp = othereast_temp + data.nuv1*psik*Flux*psil*weight;                
                    globalspeedF = max([localspeedF globalspeedF]);

                    %West Boundary terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1-data.deltav1/2;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                    quadpoint = [v1 v2 v3];

                    psik = data.vectphi_west_Trans(:,:,v1quad,v2quad);
                    psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    Flux = predictor_F_dqr(NaN,NaN,quadpoint,data.appdata);
                    fluxwest_temp = fluxwest_temp + data.nuv1*psik*Flux*psil*weight;                
                    psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
                    Flux = predictor_F_dql(NaN,NaN,quadpoint,data.appdata);
                    otherwest_temp = otherwest_temp + data.nuv1*psik*Flux*psil*weight;
                    globalspeedF = max([localspeedF globalspeedF]);

                    if data.space_dims >= 2
                        %nort Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                        quadpoint = [v1 v2 v3];

                        psik = data.vectphi_nort_Trans(:,:,v1quad,v2quad);
                        psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        Flux = predictor_G_dql(NaN,NaN,quadpoint,data.appdata);
                        fluxnort_temp = fluxnort_temp + data.nuv2*psik*Flux*psil*weight;                
                        psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                        Flux = predictor_G_dqr(NaN,NaN,quadpoint,data.appdata);
                        othernort_temp = othernort_temp + data.nuv2*psik*Flux*psil*weight;                
                        globalspeedG = max([localspeedG globalspeedG]);

                        %sout Boundary terms
                        v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                        v2 = data.v2centers(iv2)+data.deltav2*ir2-data.deltav2/2;         
                        v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                        quadpoint = [v1 v2 v3];

                        psik = data.vectphi_sout_Trans(:,:,v1quad,v2quad);
                        psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        Flux = predictor_G_dqr(NaN,NaN,quadpoint,data.appdata);
                        fluxsout_temp = fluxsout_temp + data.nuv2*psik*Flux*psil*weight;                
                        psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                        Flux = predictor_G_dql(NaN,NaN,quadpoint,data.appdata);
                        othersout_temp = othersout_temp + data.nuv2*psik*Flux*psil*weight;
                        globalspeedG = max([localspeedG globalspeedG]);
                    end
                end

                for k = 1:data.Pdp1
                    tquad = data.Dp1list(k,1); 
                    %tloc = data.Dp1quadlocs(k,1);
                    v1quad = data.Dp1list(k,2); 
                    v1loc = data.Dp1quadlocs(k,2);
                    if data.space_dims >= 2
                        v2quad = data.Dp1list(k,3);
                        v2loc = data.Dp1quadlocs(k,3);
                    else
                        v2quad = 1;
                        v2loc = inf;
                    end
                    if data.space_dims >= 3
                        v3quad = data.Dp1list(k,4);
                        v3loc = data.Dp1quadlocs(k,4);
                    else
                        v3quad = 1;
                        v3loc = inf;
                    end
                    weight = data.Dp1quadwgts(k);
                    %Derivative terms
                    v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                    v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v2loc;         
                    v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v3loc;   
                    quadpoint = [v1 v2 v3];

                    phil = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
                    psikdxii = data.vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
                    [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
                    %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                    globalspeedF = max([localspeedF globalspeedF]);
                    fterm = weight*data.nuv1*psikdxii*f*phil;
                    jaccell_temp = jaccell_temp + fterm;
                    if data.space_dims >= 2
                        phikdeta = data.vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
                        [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                        %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                        globalspeedG = max([localspeedG globalspeedG]);
                        gterm = weight*data.nuv2*phikdeta*g*phil;
                        jaccell_temp = jaccell_temp + gterm;
                    end               
                end

                RHS = fluxeast_temp - fluxwest_temp - jaccell_temp;            
                if data.space_dims >= 2
                    RHS = RHS + fluxnort_temp - fluxsout_temp;
                end           
                RHS = [ othereast_temp , RHS , -otherwest_temp ];
                if data.space_dims >= 2
                    RHS = [ othernort_temp , RHS , -othersout_temp ];
                end 
           data.corrector_update(:,:,iv1,iv2,iv3) = -RHS/2^(data.space_dims);        
        end
        end
        end
    end
end

data.Psi_futr = data.I_futr - data.Psi_dtau;
