function [ DGprediction,auxiliary,speedmaxF,speedmaxG,speedmaxH] = predictor_NewtonIteration_3D(DGpast,auxiliary,data)

speedmaxF = 0;
speedmaxG = 0;
speedmaxH = 0;

deltar = -data.r_param:data.r_param;

thetaT = data.thetaT;
Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
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
nuv3 = data.nuv3;

predictor_quadrature = data.predictor_quadrature;

for iv1 = 1:Nv1
    for iv2 = 1:Nv2
        for iv3 = 1:Nv3       
            region_iv1s = periodic(iv1 + deltar,Nv1);
            region_iv2s = periodic(iv2 + deltar,Nv2);
            region_iv3s = periodic(iv3 + deltar,Nv3);
            DGregion = DGpast( :, region_iv1s , region_iv2s, region_iv3s );
            DGregion_past = DGpast( :, region_iv1s , region_iv2s , region_iv3s  );
            iters = 0;        
            check = 1;                     
            cellindex_center = [iv1 iv2 iv3];
            while check            
                switch predictor_quadrature
                    case 'quadrature'
                        [Jacobian,residual,auxiliary] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data,cellindex_center);
                    case 'quadrature_free_volumes'
                        [Jacobian,residual,auxiliary] = problem_region_Jacobian_residual_quadrature(DGregion,DGregion_past,auxiliary,data,cellindex_center);
                    case 'quadrature_free_all'
                        error(' Quadrature routine not implemented.')                                 
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

                searchdir(:,1,1,1) = soln(ind(1));
                searchdir(:,2,1,1) = soln(ind(2));
                searchdir(:,3,1,1) = soln(ind(3));
                searchdir(:,1,2,1) = soln(ind(4));
                searchdir(:,2,2,1) = soln(ind(5));
                searchdir(:,3,2,1) = soln(ind(6));
                searchdir(:,1,3,1) = soln(ind(7));
                searchdir(:,2,3,1) = soln(ind(8));
                searchdir(:,3,3,1) = soln(ind(9));
                
                searchdir(:,1,1,2) = soln(ind(10));
                searchdir(:,2,1,2) = soln(ind(11));
                searchdir(:,3,1,2) = soln(ind(12));
                searchdir(:,1,2,2) = soln(ind(13));
                searchdir(:,2,2,2) = soln(ind(14));
                searchdir(:,3,2,2) = soln(ind(15));
                searchdir(:,1,3,2) = soln(ind(16));
                searchdir(:,2,3,2) = soln(ind(17));
                searchdir(:,3,3,2) = soln(ind(18));

                searchdir(:,1,1,3) = soln(ind(19));
                searchdir(:,2,1,3) = soln(ind(20));
                searchdir(:,3,1,3) = soln(ind(21));
                searchdir(:,1,2,3) = soln(ind(22));
                searchdir(:,2,2,3) = soln(ind(23));
                searchdir(:,3,2,3) = soln(ind(24));
                searchdir(:,1,3,3) = soln(ind(25));
                searchdir(:,2,3,3) = soln(ind(26));
                searchdir(:,3,3,3) = soln(ind(27));
              
                DGregion = DGregion - searchdir;
            end
            DGprediction(:,iv1,iv2,iv3) = DGregion(:,2,2,2);
        end
    end
end

end

