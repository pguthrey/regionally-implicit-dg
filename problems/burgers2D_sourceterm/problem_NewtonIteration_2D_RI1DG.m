function [DGprediction,speedmaxF,speedmaxG,residuals] = problem_NewtonIteration_2D_RI1DG(DGnm1,DGprev,data,residuals)

DGprediction = zeros(data.thetaT,data.Nv1,data.Nv2);
        speedmaxF = 0;
        speedmaxG = 0;
for i = 1:data.Nv1
    for j = 1:data.Nv2
        
        Ip1 = mod(i+1-1,data.Nv1)+1;
        Im1 = mod(i-1-1,data.Nv1)+1;
        
        Jp1 = mod(j+1-1,data.Nv2)+1;
        Jm1 = mod(j-1-1,data.Nv2)+1;
        
        q_sw = DGnm1(:,Im1,Jm1);
        q_sc = DGnm1(:,i,Jm1);
        q_se = DGnm1(:,Ip1,Jm1);

        q_cw = DGnm1(:,Im1,j);
        q_cc = DGnm1(:,i,j);
        q_ce = DGnm1(:,Ip1,j);

        q_nw = DGnm1(:,Im1,Jp1);
        q_nc = DGnm1(:,i,Jp1);
        q_ne = DGnm1(:,Ip1,Jp1);        
        
        q_sw_past = DGprev(:,Im1,Jm1);
        q_sc_past = DGprev(:,i,Jm1);
        q_se_past = DGprev(:,Ip1,Jm1);

        q_cw_past = DGprev(:,Im1,j);
        q_cc_past = DGprev(:,i,j);
        q_ce_past = DGprev(:,Ip1,j);

        q_nw_past = DGprev(:,Im1,Jp1);
        q_nc_past = DGprev(:,i,Jp1);
        q_ne_past = DGprev(:,Ip1,Jp1);        
        
        
        switch data.M
            case 1
                [DGprediction_cell,speedmaxF_cell,speedmaxG_cell] = problem_predictor_M1(q_sw,q_sc,q_se, ...
                                                                            q_cw,q_cc,q_ce, ...
                                                                            q_nw,q_nc,q_ne, ...
                                                                            q_sw_past,q_sc_past,q_se_past, ...
                                                                            q_cw_past,q_cc_past,q_ce_past, ...
                                                                            q_nw_past,q_nc_past,q_ne_past);        
            case 2
                [DGprediction_cell,speedmaxF_cell,speedmaxG_cell] = problem_predictor_M2(q_sw,q_sc,q_se, ...
                                                                            q_cw,q_cc,q_ce, ...
                                                                            q_nw,q_nc,q_ne, ...
                                                                            q_sw_past,q_sc_past,q_se_past, ...
                                                                            q_cw_past,q_cc_past,q_ce_past, ...
                                                                            q_nw_past,q_nc_past,q_ne_past,data);        
            case 4
                [DGprediction_cell,speedmaxF_cell,speedmaxG_cell] = problem_predictor_M4(q_sw,q_sc,q_se, ...
                                                                            q_cw,q_cc,q_ce, ...
                                                                            q_nw,q_nc,q_ne, ...
                                                                            q_sw_past,q_sc_past,q_se_past, ...
                                                                            q_cw_past,q_cc_past,q_ce_past, ...
                                                                            q_nw_past,q_nc_past,q_ne_past);        
        end
        speedmaxF = max(speedmaxF,speedmaxF_cell);
        speedmaxG = max(speedmaxG,speedmaxG_cell);
        DGprediction(:,i,j) = DGprediction_cell;
        
    end
end

