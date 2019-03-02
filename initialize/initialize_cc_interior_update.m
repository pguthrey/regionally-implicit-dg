function [update,update_index,Kmax] = initialize_cc_interior_update(data)
% written by Pierson Guthrey

index = @(gamma) (gamma-1)*data.thetaT+(1:data.thetaT);
regions = 2*data.r_param+1;
matsize = regions^data.space_dims*data.thetaT;
r = data.r_param;

predictor_update = NaN(data.thetaT,matsize,data.Nv1,data.Nv2,data.Nv3);

ry = 0;
rz = 0;
if data.space_dims >= 2
    ry = data.r_param;
    if data.space_dims >= 3
        rz = data.r_param;
    end
end

Z = [eye(data.theta);zeros(data.thetaT-data.theta,data.theta)];

Z = 0;
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
            phik = data.vectpsi_Trans(:,:,tquad,v1quad,v2quad);
            phil = data.vectphi(:,:,v1quad,v2quad);
            Z = Z + phik*phil*weight/2^(data.space_dims+1);
end

scale = @(ir1,ir2,ir3) ir1 + r + 1 + (2*r+1)*(ir2 + ry) + (2*r+1)*(2*ry+1)*(ir3 + rz)   
row = @(ir1,ir2,ir3) (1:data.thetaT)+ (scale(ir1,ir2,ir3)-1)*data.thetaT;
col = @(ir1,ir2,ir3) (1:data.thetaT)+ (scale(ir1,ir2,ir3)-1)*data.thetaT;

iv1vals = [1 floor(data.Nv1/2) data.Nv1];
for iv1_k = 1:3
    iv1_
for iv2 = 1
for iv3 = 1
    LHS = zeros(matsize);
    RHS = LHS;
    iv1 = iv1vals(iv1_k);
    for ir3 = -rz:rz            
    for ir2 = -ry:ry
    for ir1 = -r:r
      
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
            Flux = compute_F_dflux_dql(0,0,quadpoint,data.appdata);
            fluxeast_temp = fluxeast_temp + data.nuv1*psik*Flux*psil*weight;
            [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
            trunceast_temp = trunceast_temp + data.nuv1*psik*f*psil*weight;
            psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);
            Flux = compute_F_dflux_dqr(0,0,quadpoint,data.appdata);
            othereast_temp = othereast_temp + data.nuv1*psik*Flux*psil*weight;   

            %West Boundary terms
            v1 = data.v1centers(iv1)+data.deltav1*ir1-data.deltav1/2;        
            v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2*v1loc;         
            v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
            quadpoint = [v1 v2 v3];

            psik = data.vectpsi_west_Trans(:,:,tquad,v1quad,v2quad);
            psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
            Flux = compute_F_dflux_dqr(0,0,quadpoint,data.appdata);
            fluxwest_temp = fluxwest_temp + data.nuv1*psik*Flux*psil*weight;
            [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
            truncwest_temp = truncwest_temp + data.nuv1*psik*f*psil*weight;
            psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad); 
            Flux = compute_F_dflux_dql(0,0,quadpoint,data.appdata);
            otherwest_temp = otherwest_temp + data.nuv1*psik*Flux*psil*weight;
            if data.space_dims >= 2
                %nort Boundary terms
                v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                v2 = data.v2centers(iv2)+data.deltav2*ir2+data.deltav2/2;         
                v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;   
                quadpoint = [v1 v2 v3];

                psik = data.vectpsi_nort_Trans(:,:,tquad,v1quad,v2quad);
                psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dql(0,0,quadpoint,data.appdata);
                fluxnort_temp = fluxnort_temp + data.nuv2*psik*Flux*psil*weight;                
                [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                truncnort_temp = truncnort_temp + data.nuv2*psik*g*psil*weight;
                psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dqr(0,0,quadpoint,data.appdata);
                othernort_temp = othernort_temp + data.nuv2*psik*Flux*psil*weight;                

                %sout Boundary terms
                v1 = data.v1centers(iv1)+data.deltav1*ir1+data.deltav1/2*v1loc;        
                v2 = data.v2centers(iv2)+data.deltav2*ir2-data.deltav2/2;         
                v3 = data.v3centers(iv3)+data.deltav3*ir3+data.deltav3/2*v2loc;    
                quadpoint = [v1 v2 v3];

                psik = data.vectpsi_sout_Trans(:,:,tquad,v1quad,v2quad);
                psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dqr(0,0,quadpoint,data.appdata);
                fluxsout_temp = fluxsout_temp + data.nuv2*psik*Flux*psil*weight;                
                [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                truncsout_temp = truncsout_temp + data.nuv2*psik*g*psil*weight;
                psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dql(0,0,quadpoint,data.appdata);
                othersout_temp = othersout_temp + data.nuv2*psik*Flux*psil*weight;
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
            fterm = weight*data.nuv1*psikdxii*f*phil;
            jaccell_temp = jaccell_temp + fterm;
            if data.space_dims >= 2
                psikdeta = data.vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
                [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                gterm = weight*data.nuv2*psikdeta*g*phil;
                jaccell_temp = jaccell_temp + gterm;
            end               
        end
        LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = data.I_futr - data.Psi_dtau -  jaccell_temp;
        RHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = data.I_past;
        %West LHS flux
        if ir1 > -r
            LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) - fluxwest_temp;
            LHS(row(ir1,ir2,ir3),col(ir1-1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1-1,ir2,ir3)) - otherwest_temp;
        else
            LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) - truncwest_temp;
        end
        %East LHS flux
        if ir1 < r
            LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + fluxeast_temp;
            LHS(row(ir1,ir2,ir3),col(ir1+1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1+1,ir2,ir3)) + othereast_temp;
        else
            LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + trunceast_temp;
        end            
        if data.space_dims >= 2
            %sout LHS flux
            if ir2 > -ry
                LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) - fluxsout_temp;
                LHS(row(ir1,ir2,ir3),col(ir1,ir2-1,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2-1,ir3)) - othersout_temp;
            else
                LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) - truncsout_temp;
            end
            %nort LHS flux
            if ir2 < ry
                LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + fluxnort_temp;
                LHS(row(ir1,ir2,ir3),col(ir1,ir2+1,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2+1,ir3)) + othernort_temp;
            else
                LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + truncnort_temp;
            end            
            if data.space_dims >= 3
                %down LHS flux
                if ir3 > -rz
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + this.Jac_down_cell(:,:,celli,cellj,cellk);
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3-1)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3-1)) + this.Jac_down_other(:,:,celli,cellj,cellk);
                else
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + this.Jac_trunc_down(:,:,celli,cellj,cellk);
                end
                %uppr LHS flux
                if ir3 < rz
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + this.Jac_uppr_cell(:,:,celli,cellj,cellk);
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3+1)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3+1)) + this.Jac_uppr_other(:,:,celli,cellj,cellk);
                else
                    LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) = LHS(row(ir1,ir2,ir3),col(ir1,ir2,ir3)) + this.Jac_trunc_uppr(:,:,celli,cellj,cellk);
                end
            end
        end
        %iterate cell number
        %cellnum = cellnum + 1;
    end
    end
    end
            
    coeffs = index(data.main_cell);    
    if data.smartsolver
        [ prediction ] = smartsolver(LHS,eye(size(LHS)),coeffs);
    else
        solution = inv(LHS);
        prediction = solution(coeffs,:); 
    end
    predictor_update(:,:,iv1,iv2,iv3) = prediction;
end
end
end
corrector_update =  zeros(data.theta,data.thetaT,1+2*data.space_dims,data.Nv1,data.Nv2,data.Nv3);
corrector_update_index =  zeros(3,1+2*data.space_dims,data.Nv1,data.Nv2,data.Nv3);

for iv1_k = iv1vals
    iv1 = iv1vals(iv1_k);
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
            v1 = data.v1centers(iv1)+data.deltav1/2;        
            v2 = data.v2centers(iv2)+data.deltav2/2*v1loc;         
            v3 = data.v3centers(iv3)+data.deltav3/2*v2loc;   
            quadpoint = [v1 v2 v3];

            psik = data.vectphi_east_Trans(:,:,v1quad,v2quad);
            psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
            %[f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
            %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
            Flux = compute_F_dflux_dql(0,0,quadpoint,data.appdata);
            fluxeast_temp = fluxeast_temp + data.nuv1*psik*Flux*psil*weight;                
            psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
            Flux = compute_F_dflux_dqr(0,0,quadpoint,data.appdata);
            othereast_temp = othereast_temp + data.nuv1*psik*Flux*psil*weight;                

            %West Boundary terms
            v1 = data.v1centers(iv1)-data.deltav1/2;        
            v2 = data.v2centers(iv2)+data.deltav2/2*v1loc;         
            v3 = data.v3centers(iv3)+data.deltav3/2*v2loc;    
            quadpoint = [v1 v2 v3];

            psik = data.vectphi_west_Trans(:,:,v1quad,v2quad);
            psil = data.vectpsi_west(:,:,tquad,v1quad,v2quad);   
            %[f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
            %localspeedF = problem_F_Jacobianspectralradius(NaN,quadpoint,data.appdata);
            Flux = compute_F_dflux_dqr(0,0,quadpoint,data.appdata);
            fluxwest_temp = fluxwest_temp + data.nuv1*psik*Flux*psil*weight;                
            psil = data.vectpsi_east(:,:,tquad,v1quad,v2quad);   
            Flux = compute_F_dflux_dql(0,0,quadpoint,data.appdata);
            otherwest_temp = otherwest_temp + data.nuv1*psik*Flux*psil*weight;

            if data.space_dims >= 2
                %nort Boundary terms
                v1 = data.v1centers(iv1)+data.deltav1/2*v1loc;        
                v2 = data.v2centers(iv2)+data.deltav2/2;         
                v3 = data.v3centers(iv3)+data.deltav3/2*v2loc;   
                quadpoint = [v1 v2 v3];

                psik = data.vectphi_nort_Trans(:,:,v1quad,v2quad);
                psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                %[g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                Flux = compute_G_dflux_dql(0,0,quadpoint,data.appdata);
                fluxnort_temp = fluxnort_temp + data.nuv2*psik*Flux*psil*weight;                
                psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dqr(0,0,quadpoint,data.appdata);
                othernort_temp = othernort_temp + data.nuv2*psik*Flux*psil*weight;                

                %sout Boundary terms
                v1 = data.v1centers(iv1)+data.deltav1/2*v1loc;        
                v2 = data.v2centers(iv2)-data.deltav2/2;         
                v3 = data.v3centers(iv3)+data.deltav3/2*v2loc;    
                quadpoint = [v1 v2 v3];
                psik = data.vectphi_sout_Trans(:,:,v1quad,v2quad);
                psil = data.vectpsi_sout(:,:,tquad,v1quad,v2quad);   
                %[g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                %localspeedG = problem_G_Jacobianspectralradius(NaN,quadpoint,data.appdata);
                Flux = compute_G_dflux_dqr(0,0,quadpoint,data.appdata);
                fluxsout_temp = fluxsout_temp + data.nuv2*psik*Flux*psil*weight;                
                psil = data.vectpsi_nort(:,:,tquad,v1quad,v2quad);   
                Flux = compute_G_dflux_dql(0,0,quadpoint,data.appdata);
                othersout_temp = othersout_temp + data.nuv2*psik*Flux*psil*weight;
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
            v1 = data.v1centers(iv1)+data.deltav1/2*v1loc;        
            v2 = data.v2centers(iv2)+data.deltav2/2*v2loc;         
            v3 = data.v3centers(iv3)+data.deltav3/2*v3loc;   
            quadpoint = [v1 v2 v3];

            phil = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
            psikdxii = data.vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
            [f] = problem_F_FluxJacobian(NaN,quadpoint,data.appdata);
            fterm = weight*data.nuv1*psikdxii*f*phil;
            jaccell_temp = jaccell_temp + fterm;
            if data.space_dims >= 2
                phikdeta = data.vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
                [g] = problem_G_FluxJacobian(NaN,quadpoint,data.appdata);
                gterm = weight*data.nuv2*phikdeta*g*phil;
                jaccell_temp = jaccell_temp + gterm;
            end               
        end

        switch data.space_dims
            case 1
                corrector_update(:,:,1,iv1,iv2,iv3) = -otherwest_temp;
                corrector_update(:,:,2,iv1,iv2,iv3) = fluxeast_temp - fluxwest_temp - jaccell_temp;
                corrector_update(:,:,3,iv1,iv2,iv3) = othereast_temp;
                
                corrector_update_index(:,1,iv1,iv2,iv3) = [iv1-1 0 0];
                corrector_update_index(:,2,iv1,iv2,iv3) = [iv1 0 0];
                corrector_update_index(:,3,iv1,iv2,iv3) = [iv1+1 0 0];
                
            case 2
                corrector_update(:,:,1,iv1,iv2,iv3) = -othersout_temp;
                corrector_update(:,:,2,iv1,iv2,iv3) = -otherwest_temp;
                corrector_update(:,:,3,iv1,iv2,iv3) = fluxnort_temp - fluxsout_temp + fluxeast_temp - fluxwest_temp - jaccell_temp;
                corrector_update(:,:,4,iv1,iv2,iv3) = othereast_temp;
                corrector_update(:,:,5,iv1,iv2,iv3) = othernort_temp;
                
                corrector_update_index(:,1,iv1,iv2,iv3) = [iv1 iv2-1 0];
                corrector_update_index(:,2,iv1,iv2,iv3) = [iv1-1 iv2 0];
                corrector_update_index(:,3,iv1,iv2,iv3) = [iv1 iv2 0];
                corrector_update_index(:,4,iv1,iv2,iv3) = [iv1+1 iv2 0];
                corrector_update_index(:,5,iv1,iv2,iv3) = [iv1 iv2+1 0];
                
            case 3
                error('whoops')
                
                
        end
end
end
end

corrector_update = corrector_update/2^data.space_dims;

rx = data.r_param;
cbar_x = 2+rx;
cbar_y = 1;
cbar_z = 1;
deltaijk = [[-1 0 0 ];[0 0 0];[1 0 0]];
if data.space_dims >= 2
    deltaijk = [[0 -1 0]; deltaijk ; [0 1 0]];
    cbar_y = 2+ry;
    if data.space_dims >= 3
        deltaijk = [[0 0 -1]; deltaijk ; [0 0 1]];
        cbar_z = 2+rz;
    end
end

Kmax = (2*cbar_x-1)*(2*cbar_y-1)*(2*cbar_z-1);
update = zeros(data.theta,data.theta,Kmax,data.Nv1,data.Nv2,data.Nv3);
update_index = NaN(3,Kmax,data.Nv1,data.Nv2,data.Nv3);


switch data.space_dims
    case 1
    update_K = @(iux,iuy,iuz) 2+rx+iux;                
    case 2
    update_K = @(iux,iuy,iuz) 2+rx+iux ...
                        + (3+2*rx)*(2+ry+iuy-1);        
    case 3
    update_K = @(iux,iuy,iuz) 2+rx+iux ...
                        + (3+2*rx)*(2+ry+iuy-1) ...
                        + (3+2*rx)*(3+2*ry)*(2+rz+iuz-1);
end

thisy_up = 1;
thisz_up = 1;
switch data.space_dims 
    case 1
        for iv1 = 1:data.Nv1
            for iux = -(rx+1):(rx+1)           
                thisx_up = mod(iv1 + iux -1,data.Nv1)+1;
                uk = update_K(iux,NaN,NaN);
                update_index(:,uk,iv1) = [thisx_up thisy_up thisz_up];
            end
        end
    case 2
        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2
            for iuy = -(ry+1):(ry+1)
            for iux = -(rx+1):(rx+1)           
                thisx_up = mod(iv1 + iux -1,data.Nv1)+1;
                thisy_up = mod(iv2 + iuy -1,data.Nv2)+1;
                uk = update_K(iux,iuy,NaN);
                update_index(:,uk,iv1,iv2) = [thisx_up thisy_up thisz_up];
            end
            end
        end
        end
    case 3
        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2
        for iv3 = 1:data.Nv3
            for iuz = -(rz+1):(rz+1)
            for iuy = -(ry+1):(ry+1)
            for iux = -(rx+1):(rx+1)           
                thisx_up = mod(iv1 + iux -1,data.Nv1)+1;
                thisy_up = mod(iv2 + iuy -1,data.Nv2)+1;
                thisz_up = mod(iv3 + iuz -1,data.Nv3)+1;            
                uk = update_K(iux,iuy,iuz);
                update_index(:,uk,iv1,iv2,iv3) = [thisx_up thisy_up thisz_up];
            end
            end
            end
        end
        end
        end
end

Ipast = data.I_past;
for iv1_k = 1:length(iv1vals)
for iv2 = 1:data.Nv2
for iv3 = 1:data.Nv3
    for delta_ci = 1:length(deltaijk(:,1))
        C = corrector_update(:,:,delta_ci,iv1_k,iv2,iv3);
        delta_i = deltaijk(delta_ci,1);
        delta_j = deltaijk(delta_ci,2);
        delta_k = deltaijk(delta_ci,3);
        thisx = mod(iv1_k + delta_i - 1,data.Nv1) + 1;
        thisy = mod(iv2 + delta_j - 1,data.Nv2) + 1;
        thisz = mod(iv3 + delta_k - 1,data.Nv3) + 1;
        for irz = -rz:rz
        for iry = -ry:ry
        for irx = -rx:rx
            %{
            pred_K = 1+rx+irx+ ...
                     (1+2*rx)*(1+ry+iry-1) ...
                     + (1+2*rx)*(1+2*ry)*(1+rz+irz-1);                 
            P = predictor_update(:,index(pred_K),thisx,thisy,thisz);
            %}
            
            P = predictor_update(:,col(irx,iry,irz),thisx,thisy,thisz);            
            
            %{
            thisx_up = mod(thisx + irx -1,data.Nv1)+1;
            thisy_up = mod(thisy + iry -1,data.Nv2)+1;
            thisz_up = mod(thisz + irz -1,data.Nv3)+1;
            %}
            
            iux = irx + delta_i;
            iuy = iry + delta_j;
            iuz = irz + delta_k;
            uk = update_K(iux,iuy,iuz);
            update(:,:,uk,iv1,iv2,iv3) = update(:,:,uk,iv1,iv2,iv3) + C*P*Ipast*Z;

        end
        end
        end
    end
end
end
end

end

