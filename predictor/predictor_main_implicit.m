function [ DGprediction,auxiliary,speedmaxF,speedmaxG,speedmaxH] = predictor_main_implicit(DGpast_phi,auxiliary,data)
% written by Pierson Guthrey


[ DGpast ] = predictor_project_phi2psi(DGpast_phi,data);

%{
DGinit = DGpast;
Nv1 = data.Nv1;
Nv2 = data.Nv2;
nuv1 = data.nuv1;
nuv2 = data.nuv2;

for iv1 = 1:Nv1
    for iv2 = 1:Nv2
        DGcell = DGpast(:,iv1,iv2);
        qstar = DGcell;
        qstar_past = DGcell;
        switch data.predictorbasis
            case 'P'
                switch data.M
                    case 4
                        [Jac_cell,residual_cell] = problem_exact_jaccell_M4_P(qstar,qstar_past,nuv1,nuv2);
                    case 6
                        [Jac_cell,residual_cell] = problem_exact_jaccell_M6_P(qstar,qstar_past,nuv1,nuv2);
                    case 8
                        [Jac_cell,residual_cell] = problem_exact_jaccell_M8_P(qstar,qstar_past,nuv1,nuv2);
                    otherwise
                        error('This case is not yet implemented')
                end           
            otherwise
                error('This case is not yet implemented')
        end                      
        [east_trunc,west_trunc, ...
        Jac_east_trunc,Jac_west_trunc, ...
        nort_trunc,sout_trunc, ...
        Jac_nort_trunc,Jac_sout_trunc, ...
        maxspeedF,maxspeedG] ...
        = problem_precompute_local_all(DGcell,data,1);
        
        residual = residual_cell + east_trunc + west_trunc ...
                                 + nort_trunc + sout_trunc; 
    
        Jacobian = Jac_cell + Jac_east_trunc + Jac_west_trunc ...
                            + Jac_nort_trunc + Jac_sout_trunc; 
        DGinit(:,iv1,iv2) = DGcell - Jacobian\residual;
    end
end

DGpast = DGinit;
%}

speedmaxG = 0;
speedmaxH = 0;        
switch data.space_dims
    case 1
        [DGprediction,auxiliary,speedmaxF]                     = predictor_NewtonIteration_1D(DGpast,auxiliary,data);
    case 2
        [DGprediction,auxiliary,speedmaxF,speedmaxG]           = predictor_NewtonIteration_2D(DGpast,auxiliary,data);
    case 3
        [DGprediction,auxiliary,speedmaxF,speedmaxG,speedmaxH] = predictor_NewtonIteration_3D(DGpast,auxiliary,data);
otherwise
    error('invalid space dimensions')
end
    

%{       
close all

subplot(1,3,1)
A(:,:) = max(abs(DGprediction-DGpast),[],1)
surf(A)
colorbar

subplot(1,3,2)
A(:,:) = max(abs(DGinit-DGpast),[],1)
surf(A)
colorbar


subplot(1,3,3)
A(:,:) = max(abs(DGinit-DGprediction),[],1)
surf(A)
colorbar

keyboard     

%}
end
