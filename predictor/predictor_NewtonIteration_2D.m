function [ DGprediction,auxiliary,speedmaxF,speedmaxG] = predictor_NewtonIteration_2D(DGpast,auxiliary,data)



speedmaxF = 0;
speedmaxG = 0;

deltar = -data.r_param:data.r_param;

thetaT = data.thetaT;
Nv1 = data.Nv1;
Nv2 = data.Nv2;
%{
solverinfo.thetaT = data.thetaT;
solverinfo.cells_per_region = data.cells_per_region;
solverinfo.region_per_dimension = data.region_per_dimension;
solverinfo.space_dims = data.space_dims;
solverinfo.main_cell = data.main_cell;
solverinfo.rx_param = data.rx_param;
solverinfo.ry_param = data.ry_param;
solverinfo.rz_param = data.rz_param;
solverinfo.smartsolver = data.smartsolver;
%}

periodic = @(i,N) mod(i-1,N)+1;
DGprediction = DGpast;
ind = @(i) (1:data.thetaT) + data.thetaT*(i-1);

DGregion = NaN(thetaT,3,3);
searchdir = NaN(thetaT,3,3);

nuv1 = data.nuv1;
nuv2 = data.nuv2;

predictor_quadrature = data.predictor_quadrature;

for iv1 = 1:Nv1
    for iv2 = 1:Nv2
       
        region_iv1s = periodic(iv1 + deltar,Nv1);
        region_iv2s = periodic(iv2 + deltar,Nv2);
        DGregion(:,:,:) = DGpast( :, region_iv1s , region_iv2s );
        DGregion_past(:,:,:) = DGpast( :, region_iv1s , region_iv2s );
        iters = 0;        
        check = 1;                      
        while check            
            switch predictor_quadrature
                case 'quadrature'
                    [Jacobian,residual] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data);
                case 'quadrature_free_volumes'
                    [Jacobian,residual] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data);
                case 'quadrature_free_all'
                    switch data.M
                        case 2
                            switch data.predictorbasis
                                case 'P'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M2_P(DGregion,DGregion_past,nuv1,nuv2); 
                                case 'Q'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M2_Q(DGregion,DGregion_past,nuv1,nuv2); 
                                otherwise 
                                error(' Quadrature routine not implemented.')                                 
                            end
                        case 4
                            switch data.predictorbasis
                                case 'P'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M4_P(DGregion,DGregion_past,nuv1,nuv2); 
                                case 'Q'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M4_Q(DGregion,DGregion_past,nuv1,nuv2); 
                                otherwise 
                                error(' Quadrature routine not implemented.')                                 
                            end
                        case 6
                            switch data.predictorbasis
                                case 'P'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M6_P(DGregion,DGregion_past,nuv1,nuv2); 
                                case 'Q'
                                [Jacobian,residual] = problem_region_Jacobian_residual_M6_Q(DGregion,DGregion_past,nuv1,nuv2); 
                            end
                        otherwise
                            error(' Quadrature routine not implemented.')
                    end                
                otherwise
                    error(' Quadrature routine not implemented.')
            end

            soln = Jacobian\residual;

            this_residual_vec = residual(4*data.thetaT+(1:data.thetaT));
            this_residual = max(abs(this_residual_vec));
                        
            iters = iters + 1;
            
            check1 = this_residual > data.residTol; 
            check2 = iters < data.maxiters;
            check = check1 && check2 ;
%{
            [Jacobian_test,residual_test] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,data); 
            soln_test = Jacobian_test\residual_test;
            Jacerror = max(max(abs(Jacobian_test - Jacobian)));
            reserror = max(abs(residual_test - residual));
            solerror = max(abs(soln_test - soln));
            
            if max([Jacerror reserror solerror]) > 1e-6
                Jacerror = Jacerror 
                reserror = reserror 
                solerror = solerror 
                
                figure(1)
                subplot(1,3,1)
                surf(Jacobian_test)
                shading flat
                view(0,-90)
                colorbar
                subplot(1,3,2)
                surf(Jacobian)
                shading flat
                view(0,-90)
                colorbar
                subplot(1,3,3)
                surf(log10(abs(Jacobian_test-Jacobian)+2*eps))
                shading flat
                axis([1 data.thetaT*9 1 data.thetaT*9 -5 1 -5 1])
                view(0,-90)
                colorbar

                ymin = min([residual_test;residual]);
                ymax = max([residual_test;residual]);

                figure(2)
                subplot(1,3,1)
                plot(residual_test)
                axis([1 data.thetaT*9 ymin-.05 ymax+.05])
                subplot(1,3,2)
                plot(residual) 
                axis([1 data.thetaT*9 ymin-.05 ymax+.05])
                subplot(1,3,3)
                plot(abs(residual_test-residual))


                figure(3)
                plot(diag(Jacobian_test)-diag(Jacobian))                


               keyboard  
            end
            %{



            keyboard             
            %}
            %}
            
            searchdir(:,1,1) = soln(ind(1));
            searchdir(:,2,1) = soln(ind(2));
            searchdir(:,3,1) = soln(ind(3));
            searchdir(:,1,2) = soln(ind(4));
            searchdir(:,2,2) = soln(ind(5));
            searchdir(:,3,2) = soln(ind(6));
            searchdir(:,1,3) = soln(ind(7));
            searchdir(:,2,3) = soln(ind(8));
            searchdir(:,3,3) = soln(ind(9));

            DGregion = DGregion - searchdir;
        end
        DGprediction(:,iv1,iv2) = DGregion(:,2,2);
        iters_map(iv1,iv2) = iters;
    end
end

%surf(iters_map)
%keyboard 


end

