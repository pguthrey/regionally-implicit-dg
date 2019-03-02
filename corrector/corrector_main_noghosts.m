function [DGcorrection,auxiliary,speedmaxF,speedmaxG,speedmaxH,data] = corrector_main_noghosts(DGprediction,DGpast,auxiliary,data)
% data.Perform the RIDG correction step
% written by data.Pierson Guthrey
% -------------------------------------------------
% INPUTS    DGprediction : the output from the prediction step
% DGpast       : coefficients from the previous timestep
% OUTPUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

%h1 = waitbar(0,'Correcting...','OuterPosition', [100 100 300 75]);

switch data.space_dims
    case 1
    [DGeast,DGwest] = problem_boundaryconditions_psi(DGprediction,data);
    case 2
    [DGeast,DGwest,DGnort,DGsout] = problem_boundaryconditions_psi(DGprediction,data);
    case 3
    [DGeast,DGwest,DGnort,DGsout,DGuppr,DGdown] = problem_boundaryconditions_psi(DGprediction,data);
end

v1centers = data.v1centers;
v2centers = data.v2centers; 
v3centers = data.v3centers;

correct_cell = NaN(data.theta,data.Nv1,data.Nv2,data.Nv3); 

speedmaxF = 0;
speedmaxG = 0;
speedmaxH = 0;

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;


PHI_norm_inv = data.PHI_norm_inv;

if strcmp(data.linearcase,'linear')
    for iv1 = 1:Nv1  
    for iv3 = 1:Nv3
    for iv2 = 1:Nv2
        cellnum = 1;
        qstar = DGprediction(:,iv1,iv2,iv3);
        qeast = DGeast(:,iv1,iv2,iv3);
        qwest = DGwest(:,iv1,iv2,iv3);
        Q = [qeast;qstar;qwest];
        if data.space_dims >= 2
            qnort = DGnort(:,iv1,iv2,iv3);
            qsout = DGsout(:,iv1,iv2,iv3);
            Q = [qnort;Q;qsout];            
            if data.space_dims >= 3
                quppr = DGuppr(:,iv1,iv2,iv3);
                qdown = DGdown(:,iv1,iv2,iv3);
            else
                quppr = 0.*qstar;
                qdown = 0.*qstar;
            end
        else 
            qnort = 0.*qstar;
            qsout = 0.*qstar;
        end
       correct_cell(:,iv1,iv2,iv3) = data.corrector_update(:,:,iv1,iv2,iv3)*Q;
    end
    end
    end
    
    
else
% -----------------------------------------------------
% Nonlinear Case
% -----------------------------------------------------

    periodic = @(i,N) mod(i-1,N)+1;
    
    for iv1 = 1:Nv1  
    for iv3 = 1:Nv3
    for iv2 = 1:Nv2

        qstar = DGprediction(:,iv1,iv2,iv3);
        qeast = DGprediction(:,periodic(iv1+1,Nv1),iv2,iv3);
        qwest = DGprediction(:,periodic(iv1-1,Nv1),iv2,iv3);
                             
        if data.space_dims >= 2
                        
            qstar = DGprediction(:,iv1,iv2,iv3);
            qnort = DGprediction(:,iv1,periodic(iv2+1,Nv2),iv3);
            qsout = DGprediction(:,iv1,periodic(iv2-1,Nv2),iv3);

            if data.space_dims >= 3
                quppr = DGuppr(:,iv1,iv2,iv3);
                qdown = DGdown(:,iv1,iv2,iv3);
            else
                quppr = 0.*qstar;
                qdown = 0.*qstar;
            end
        else 
            qnort = 0.*qstar;
            qsout = 0.*qstar;
        end
        
        v1center = v1centers(iv1); 
        v2center = v2centers(iv2); 
        v3center = v3centers(iv3); 
        cellcenter = [v1center v2center v3center];
        cellindex = [iv1 iv2 iv3];
        
        [correct_cell_temp,auxiliary,maxspeedF,data] = problem_correct_eastwest(qstar,qeast,qwest,auxiliary,data,cellcenter,cellindex);   
        speedmaxF = max([speedmaxF maxspeedF]);
        [correct_cell_temp,auxiliary,data] = problem_correct_source(correct_cell_temp,qstar,auxiliary,data,cellcenter,cellindex);
        if data.space_dims >= 2 
            [correct_cell_temp,auxiliary,maxspeedG,data] = problem_correct_nortsout(correct_cell_temp,qstar,qnort,qsout,auxiliary,data,cellcenter,cellindex);   
            speedmaxG = max([speedmaxG maxspeedG]);
            if data.space_dims >= 3 
                [correct_cell_temp,auxiliary,maxspeedH,data] = problem_correct_upprdown(correct_cell_temp,qstar,quppr,qdown,auxiliary,data,cellcenter,cellindex);  
                speedmaxH = max([speedmaxH maxspeedH]);
            end
        end
        
        %{ 
        [correct_cell_temp,maxspeedF,data] = correct_eastwest(qstar,qeast,qwest,data,cellcenter,iv1);   
        speedmaxF = max([speedmaxF maxspeedF]);
        if data.space_dims >= 2 
            [correct_cell_temp,maxspeedG] = correct_nortsout(correct_cell_temp,qstar,qnort,qsout,data,cellcenter);   
            speedmaxG = max([speedmaxG maxspeedG]);
            if data.space_dims >= 3 
                [correct_cell_temp,maxspeedH] = correct_upprdown(correct_cell_temp,qstar,quppr,qdown,data,cellcenter);  
                speedmaxH = max([speedmaxH maxspeedH]);
            end
        end
%}
        correct_cell(:,iv1,iv2,iv3) = PHI_norm_inv*correct_cell_temp;

    end
    end
        %waitbar(iv1/Nv1,h1,'Correcting...')
    end
end

%done = iv2/(Nv2+1);
%waitbar(done,h1)
DGcorrection = DGpast - correct_cell;%(:,1:Nv1,1:Nv2,1:Nv3);
%DGcorrection_out = DGcorrection;
%delete(h1) 
