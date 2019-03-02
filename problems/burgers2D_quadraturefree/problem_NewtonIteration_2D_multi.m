function [ DGprediction,speedmaxF,speedmaxG,residuals] = predictor_NewtonIteration_2D_multi(DGprev,DGpast,residuals,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%h1 = waitbar(0,'Predicting...','OuterPosition', [100 100 300 75]);

%See predictor_main_parfor

%DGprev = problem_boundaryconditions_psi(DGprevious_psi,data);

Nv1 = data.Nv1;
Nv2 = data.Nv2;
thetaT = data.thetaT;

speedmaxF = 0;
speedmaxG = 0;

[DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data);

deltar = -data.r_param:data.r_param;
%{
east_trunc = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
west_trunc = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
nort_trunc = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
sout_trunc = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
east_flux = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
west_flux = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
nort_flux = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
sout_flux = NaN(data.thetaT,data.Nv1,1+2*data.r_param);
residual_cell = NaN(data.thetaT,data.Nv1,1+2*data.r_param);

Jac_east_cell = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_west_cell = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_nort_cell = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_sout_cell = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_east_other = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_west_other = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_nort_other = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_sout_other = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_east_trunc = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_west_trunc = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_nort_trunc = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_sout_trunc = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);
Jac_cell = NaN(data.thetaT,data.thetaT,data.Nv1,1+2*data.r_param);

%}

if data.r_param > 0 
    for r = 1+deltar(1:end-1)  % [0 1]
        %We want to store these as a Nv1*3 array, otherwise memory is an
        %issue
        iv2prime = 1+mod(r-1,Nv2); % [Nv2 1]
        for iv1 = 1:Nv1    
                qstar = DGprev(:,iv1,iv2prime);
                qeast = DGprev_east(:,iv1,iv2prime);
                qwest = DGprev_west(:,iv1,iv2prime);
                qnort = DGprev_nort(:,iv1,iv2prime);
                qsout = DGprev_sout(:,iv1,iv2prime);
                qpast = DGpast(:,iv1,iv2prime);

                v1center = data.deltav1*(iv1-1/2)+data.v1_lb;%data.v1centers(iv1); 
                v2center = data.deltav2*(r-1/2)+data.v2_lb;%data.v2centers(iv2); 
                cellcenter = [v1center v2center 0];

                [east_trunc_temp, ...
                west_trunc_temp, ...
                 east_flux_temp, ...
                 west_flux_temp, ...
                residual_temp,maxspeedF_res] ...
                    = predictor_precompute_eastwest_residual(qstar,qpast,qeast,qwest,data,cellcenter);       

                [Jac_east_cell_temp, ...
                Jac_west_cell_temp, ...
                Jac_east_other_temp, ...
                Jac_west_other_temp, ...
                Jac_east_trunc_temp, ...
                Jac_west_trunc_temp, ...
                Jac_cell_temp,maxspeedF_jac] ...
                = predictor_precompute_eastwest_Jacobian(qstar,qeast,qwest,data,cellcenter);
                                                          
                [nort_trunc_temp, ...
                sout_trunc_temp, ...
                nort_flux_temp, ...
                sout_flux_temp, ...
                residual_cell_temp, ...
                maxspeedG_res] ...
                = predictor_precompute_nortsout_residual(qstar,qnort,qsout,residual_temp,data,cellcenter);

                [Jac_nort_cell_temp, ...
                Jac_sout_cell_temp, ...
                Jac_nort_other_temp, ...
                Jac_sout_other_temp, ...
                Jac_nort_trunc_temp, ...
                Jac_sout_trunc_temp, ...
                Jac_cell_temp, ...
                maxspeedG_jac] ...
                = predictor_precompute_nortsout_Jacobian(qstar,qnort,qsout,Jac_cell_temp,data,cellcenter);
                
            
                %Compute residual and Jacobian features
               [trunceast(:,iv1,r+data.r_param),truncwest(:,iv1,r+data.r_param), ...
                  fluxeast(:,iv1,r+data.r_param), fluxwest(:,iv1,r+data.r_param), ...
                  residual_cell(:,iv1,r+data.r_param), ...
                  Jac_east_cell(:,:,iv1,r+data.r_param),Jac_west_cell(:,:,iv1,r+data.r_param), ...
                Jac_east_other(:,:,iv1,r+data.r_param),Jac_west_other(:,:,iv1,r+data.r_param), ...
                Jac_trunc_east(:,:,iv1,r+data.r_param),Jac_trunc_west(:,:,iv1,r+data.r_param), ...
                Jac_cell(:,:,iv1,r+data.r_param),maxspeedF] = problem_precompute_eastwest_all(qstar,qeast,qwest,qpast,data,NaN);                
            
                [truncnort(:,iv1,r+data.r_param),truncsout(:,iv1,r+data.r_param), ...
                  fluxnort(:,iv1,r+data.r_param), fluxsout(:,iv1,r+data.r_param), ...
                  residual_cell(:,iv1,r+data.r_param), ...
                  Jac_nort_cell(:,:,iv1,r+data.r_param),Jac_sout_cell(:,:,iv1,r+data.r_param), ...
                Jac_nort_other(:,:,iv1,r+data.r_param),Jac_sout_other(:,:,iv1,r+data.r_param), ...
                Jac_trunc_nort(:,:,iv1,r+data.r_param),Jac_trunc_sout(:,:,iv1,r+data.r_param), ...
                Jac_cell(:,:,iv1,r+data.r_param),maxspeedF] = problem_precompute_nortsout_all(qstar,qnort,qsout,residual_cell(:,iv1,r+data.r_param),Jac_cell(:,:,iv1,r+data.r_param),data,NaN);                 
            
                speedmaxF = max([speedmaxF maxspeedF_jac maxspeedF_res]);
                speedmaxG = max([speedmaxG maxspeedG_jac maxspeedG_res]);

        end
    end

    DGprediction = NaN(thetaT,data.Nv1,data.Nv2);
    solverdata.thetaT = data.thetaT;
    solverdata.cells_per_region = data.cells_per_region;
    solverdata.region_per_dimension = data.region_per_dimension;
    solverdata.space_dims = data.space_dims;
    solverdata.main_cell = data.main_cell;
    solverdata.rx_param = data.rx_param;
    solverdata.ry_param = data.ry_param;
    solverdata.rz_param = data.rz_param;
    solverdata.smartsolver = data.smartsolver;

    for iv2 = 1:Nv2
        iv2prime = 1+mod(iv2+1-1,Nv2); %new cell [2:Nv2 1]
        for iv1 = 1:(Nv1)
                qstar = DGprev(:,iv1,iv2prime);
                qeast = DGprev_east(:,iv1,iv2prime);
                qwest = DGprev_west(:,iv1,iv2prime);
                qnort = DGprev_nort(:,iv1,iv2prime);
                qsout = DGprev_sout(:,iv1,iv2prime);
                qpast = DGpast(:,iv1,iv2prime);

                v1center = data.deltav1*(iv1-1/2)+data.v1_lb;%data.v1centers(iv1); 
                v2center = data.deltav2*(iv2-1/2)+data.v2_lb;% data.v2centers(iv2prime); 
                cellcenter = [v1center v2center 0]; 

                % 2*data.r_param+1 = 3
                [east_trunc(:,iv1,2*data.r_param+1), ...
                    west_trunc(:,iv1,2*data.r_param+1), ...
                 east_flux(:,iv1,2*data.r_param+1), ...
                 west_flux(:,iv1,2*data.r_param+1), ...
                residual_cell(:,iv1,2*data.r_param+1), ...
                maxspeedF_res] ...
                    = predictor_precompute_eastwest_residual(qstar,qpast,qeast,qwest,data,cellcenter);       

                [Jac_east_cell(:,:,iv1,2*data.r_param+1), ...
                 Jac_west_cell(:,:,iv1,2*data.r_param+1), ...
                Jac_east_other(:,:,iv1,2*data.r_param+1), ...
                Jac_west_other(:,:,iv1,2*data.r_param+1), ...
                Jac_east_trunc(:,:,iv1,2*data.r_param+1), ...
                Jac_west_trunc(:,:,iv1,2*data.r_param+1), ...
                Jac_cell(:,:,iv1,2*data.r_param+1), ...
                maxspeedF_jac] ...
                = predictor_precompute_eastwest_Jacobian(qstar,qeast,qwest,data,cellcenter);

                [nort_trunc(:,iv1,2*data.r_param+1), ...
                    sout_trunc(:,iv1,2*data.r_param+1), ...
                nort_flux(:,iv1,2*data.r_param+1), ...
                sout_flux(:,iv1,2*data.r_param+1), ...
                residual_cell(:,iv1,2*data.r_param+1), ...
                maxspeedG_res] ...
                = predictor_precompute_nortsout_residual(qstar,qnort,qsout,residual_cell(:,iv1,2*data.r_param+1),data,cellcenter);

                [Jac_nort_cell(:,:,iv1,2*data.r_param+1),  ...
                    Jac_sout_cell(:,:,iv1,2*data.r_param+1), ...
                Jac_nort_other(:,:,iv1,2*data.r_param+1), ...
                Jac_sout_other(:,:,iv1,2*data.r_param+1), ...
                Jac_nort_trunc(:,:,iv1,2*data.r_param+1), ...
                Jac_sout_trunc(:,:,iv1,2*data.r_param+1), ...
                Jac_cell(:,:,iv1,2*data.r_param+1),  ...
                maxspeedG_jac] ...
                = predictor_precompute_nortsout_Jacobian(qstar,qnort,qsout,Jac_cell(:,:,iv1,2*data.r_param+1),data,cellcenter);

                %Compute residual and Jacobian features
               [trunceast(:,iv1,2*data.r_param+1),truncwest(:,iv1,2*data.r_param+1), ...
                  fluxeast(:,iv1,2*data.r_param+1), fluxwest(:,iv1,2*data.r_param+1), ...
                  residual_cell(:,iv1,2*data.r_param+1), ...
                  Jac_east_cell(:,:,iv1,2*data.r_param+1),Jac_west_cell(:,:,iv1,2*data.r_param+1), ...
                Jac_east_other(:,:,iv1,2*data.r_param+1),Jac_west_other(:,:,iv1,2*data.r_param+1), ...
                Jac_trunc_east(:,:,iv1,2*data.r_param+1),Jac_trunc_west(:,:,iv1,2*data.r_param+1), ...
                Jac_cell(:,:,iv1,2*data.r_param+1),maxspeedF] = problem_precompute_eastwest_all(qstar,qeast,qwest,qpast,data,NaN);                
            
                [truncnort(:,iv1,2*data.r_param+1),truncsout(:,iv1,2*data.r_param+1), ...
                  fluxnort(:,iv1,2*data.r_param+1), fluxsout(:,iv1,2*data.r_param+1), ...
                  residual_cell(:,iv1,2*data.r_param+1), ...
                  Jac_nort_cell(:,:,iv1,2*data.r_param+1),Jac_sout_cell(:,:,iv1,2*data.r_param+1), ...
                Jac_nort_other(:,:,iv1,2*data.r_param+1),Jac_sout_other(:,:,iv1,2*data.r_param+1), ...
                Jac_trunc_nort(:,:,iv1,2*data.r_param+1),Jac_trunc_sout(:,:,iv1,2*data.r_param+1), ...
                Jac_cell(:,:,iv1,2*data.r_param+1),maxspeedF] = problem_precompute_nortsout_all(qstar,qnort,qsout,residual_cell(:,iv1,2*data.r_param+1),Jac_cell(:,:,iv1,2*data.r_param+1),data,NaN);                 
            
            
            
            
            
            
                speedmaxF = max([speedmaxF maxspeedF_jac maxspeedF_res]);
                speedmaxG = max([speedmaxG maxspeedG_jac maxspeedG_res]);

        end
        for iv1 = 1:(Nv1)
            qpastsolve = DGpast(:,iv1,iv2);   
            cellsi = 1+mod(iv1+deltar-1,Nv1);

            this.residual_cell = residual_cell(:,cellsi,:);
            this.west_trunc = west_trunc(:,cellsi,:);
            this.east_trunc = east_trunc(:,cellsi,:);
            this.east_flux = east_flux(:,cellsi,:);
            this.west_flux = west_flux(:,cellsi,:);
            this.Jac_cell = Jac_cell(:,:,cellsi,:);    
            this.Jac_east_cell = Jac_east_cell(:,:,cellsi,:); 
            this.Jac_west_cell = Jac_west_cell(:,:,cellsi,:); 
            this.Jac_east_other = Jac_east_other(:,:,cellsi,:); 
            this.Jac_west_other = Jac_west_other(:,:,cellsi,:); 
            this.Jac_east_trunc = Jac_east_trunc(:,:,cellsi,:);
            this.Jac_west_trunc = Jac_west_trunc(:,:,cellsi,:);
            this.nort_trunc = nort_trunc(:,cellsi,:);
            this.sout_trunc = sout_trunc(:,cellsi,:);      
            this.nort_flux = nort_flux(:,cellsi,:);
            this.sout_flux = sout_flux(:,cellsi,:);
            this.Jac_nort_cell = Jac_nort_cell(:,:,cellsi,:); 
            this.Jac_sout_cell = Jac_sout_cell(:,:,cellsi,:); 
            this.Jac_nort_other = Jac_nort_other(:,:,cellsi,:);
            this.Jac_sout_other = Jac_sout_other(:,:,cellsi,:);
            this.Jac_nort_trunc = Jac_nort_trunc(:,:,cellsi,:);
            this.Jac_sout_trunc = Jac_sout_trunc(:,:,cellsi,:);
            this.smartsolver = data.smartsolver;
            [DGprediction(:,iv1,iv2),residuals(iv1,iv2prime)] = predictor_solver(qpastsolve,this,solverdata); 
            
            
        end
                
        %Now we can shift memory around and move on to the next cell
        for iv1 = 1:Nv1            
            east_trunc(:,iv1,1:end-1) = east_trunc(:,iv1,2:end);
            west_trunc(:,iv1,1:end-1) = west_trunc(:,iv1,2:end);
            east_flux(:,iv1,1:end-1) = east_flux(:,iv1,2:end);
            west_flux(:,iv1,1:end-1) = west_flux(:,iv1,2:end);
            residual_cell(:,iv1,1:end-1) = residual_cell(:,iv1,2:end);
            Jac_east_cell(:,:,iv1,1:end-1) = Jac_east_cell(:,:,iv1,2:end);
            Jac_west_cell(:,:,iv1,1:end-1) = Jac_west_cell(:,:,iv1,2:end);
            Jac_east_other(:,:,iv1,1:end-1) = Jac_east_other(:,:,iv1,2:end);
            Jac_west_other(:,:,iv1,1:end-1) = Jac_west_other(:,:,iv1,2:end);
            Jac_east_trunc(:,:,iv1,1:end-1) = Jac_east_trunc(:,:,iv1,2:end);
            Jac_west_trunc(:,:,iv1,1:end-1) = Jac_west_trunc(:,:,iv1,2:end);
            Jac_cell(:,:,iv1,1:end-1) = Jac_cell(:,:,iv1,2:end);
            nort_trunc(:,iv1,1:end-1) = nort_trunc(:,iv1,2:end);
            sout_trunc(:,iv1,1:end-1) = sout_trunc(:,iv1,2:end);
            nort_flux(:,iv1,1:end-1) = nort_flux(:,iv1,2:end);
            sout_flux(:,iv1,1:end-1) = sout_flux(:,iv1,2:end);
            Jac_nort_cell(:,:,iv1,1:end-1) = Jac_nort_cell(:,:,iv1,2:end);
            Jac_sout_cell(:,:,iv1,1:end-1) = Jac_sout_cell(:,:,iv1,2:end);
            Jac_nort_other(:,:,iv1,1:end-1) = Jac_nort_other(:,:,iv1,2:end);
            Jac_sout_other(:,:,iv1,1:end-1) = Jac_sout_other(:,:,iv1,2:end);
            Jac_nort_trunc(:,:,iv1,1:end-1) = Jac_nort_trunc(:,:,iv1,2:end);
            Jac_sout_trunc(:,:,iv1,1:end-1) = Jac_sout_trunc(:,:,iv1,2:end);
        end    
        %waitbar(iv2/(Nv2+1+data.r_param),h1,'Predicting...')
    end
    
end

%delete(h1)
end

