function [ DGprediction,auxiliary,speedmaxF] = predictor_NewtonIteration_1D(DGpast,auxiliary,data)

speedmaxF = 0;

deltar = -data.r_param:data.r_param;

thetaT = data.thetaT;
Nv1 = data.Nv1;
solverinfo.thetaT = data.thetaT;
solverinfo.cells_per_region = data.cells_per_region;
solverinfo.region_per_dimension = data.region_per_dimension;
solverinfo.space_dims = data.space_dims;
solverinfo.main_cell = data.main_cell;
solverinfo.rx_param = data.rx_param;
solverinfo.ry_param = data.ry_param;
solverinfo.rz_param = data.rz_param;
solverinfo.smartsolver = data.smartsolver;

predictor_quadrature = data.predictor_quadrature;

periodic = @(i,N) mod(i-1,N)+1;
DGprediction = DGpast;
ind = @(i) (1:data.thetaT) + data.thetaT*(i-1);

DGregion = NaN(thetaT,3);
searchdir = NaN(thetaT,3);

nuv1 = data.nuv1;

for iv1 = 1:Nv1
    region_iv1s = periodic(iv1 + deltar,Nv1);
    DGregion(:,:) = DGpast( :, region_iv1s );
    DGregion_past(:,:) = DGpast( :, region_iv1s );
    iters = 0;        
    check = 1;        
    while check                      
        
        switch predictor_quadrature
            case 'quadrature'
                [Jacobian,residual] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data);
            case 'quadrature_free_all'
                switch data.M
                    case 4
                        switch data.predictorbasis
                            case 'P'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M4_P(DGregion,DGregion_past,auxiliary,nuv1); 
                            case 'Q'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M4_Q(DGregion,DGregion_past,auxiliary,nuv1); 
                        end
                    case 6
                        switch data.predictorbasis
                            case 'P'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M6_P(DGregion,DGregion_past,auxiliary,nuv1); 
                            case 'Q'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M6_Q(DGregion,DGregion_past,auxiliary,nuv1); 
                        end
                    case 8
                        switch data.predictorbasis
                            case 'P'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M8_P(DGregion,DGregion_past,auxiliary,nuv1); 
                            case 'Q'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M8_Q(DGregion,DGregion_past,auxiliary,nuv1); 
                        end
                    case 10
                        switch data.predictorbasis
                            case 'P'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M10_P(DGregion,DGregion_past,auxiliary,nuv1); 
                            case 'Q'
                            [Jacobian,residual] = problem_region_Jacobian_residual_M10_Q(DGregion,DGregion_past,auxiliary,nuv1); 
                        end
                    otherwise
                        error(' Quadrature routine not implemented.')
                end                
            otherwise
                error(' Quadrature routine not implemented.')
        end
        
        this_residual_vec = residual(data.thetaT+(1:data.thetaT));
        this_residual = max(abs(this_residual_vec));

        iters = iters + 1;

        check1 = this_residual > data.residTol; 
        check2 = iters < data.maxiters;
        check = check1 && check2 ;
        
        soln = Jacobian\residual;
        
        searchdir(:,1) = soln(ind(1));
        searchdir(:,2) = soln(ind(2));
        searchdir(:,3) = soln(ind(3));            

        DGregion = DGregion - searchdir;
    end
    DGprediction(:,iv1) = DGregion(:,2);        
end

end

