function [DFdn,dpdn_excess_cell] = compute_correlation(DGsolution,data)

GaussHermite_locs = data.GH_locs;
GaussHermite_wgts = data.GH_wgts1D;

GaussLegendre_locs = data.locs;
GaussLegendre_wgts = data.wgts1D;

deltav1 = data.deltav1;
deltav2 = data.deltav2;
deltav3 = data.deltav3;

[data] = initialize_gauss_hermite_quadrature(data,data.M);

DGaverages = DGsolution;
DGaverages(2:end,:,:,:) = 0;
DGfluctuations = DGsolution - DGaverages;

vectphi_Trans = data.vectphi_Trans;

DFdn = zeros(data.Ls,data.Nv1,data.Nv2,data.Nv3);
dpdn_excess_cell = zeros(data.Ls,data.Nv1,data.Nv2,data.Nv3);

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

            
%h = waitbar(done,'Computing Fcorr')
for iv1 = 1:Nv1 
    for iv2 = 1:Nv2 
        for iv3 = 1:Nv3
 
            n0_cell = DGaverages(1,iv1,iv2,iv3);
            [params] = cell_average_to_params(n0_cell,data);
            [x_vals,c_vals] = HNC_solve(params.Gamma,params.kappa);
            [coeff,scaling] = HNC_fit(x_vals,c_vals);


            v1bar = data.v1centers(iv1);
            v2bar = data.v1centers(iv2);
            v3bar = data.v1centers(iv3);

            quadlocsx = v1bar + data.deltav1/2*GaussLegendre_locs;
            quadlocsy = NaN;
            quadlocsz = NaN;
            if data.space_dims >= 2
                quadlocsy = v2bar + data.deltav2/2*GaussLegendre_locs;
                if data.space_dims >= 3
                    quadlocsz = v3bar + data.deltav3/2*GaussLegendre_locs;
                end
            end
                        
            DFDN_cell = 0;
            
            %{
            n0_cell = DGaverages(1,iv1,iv2,iv3);
            [params] = cell_average_to_params(n0_cell,data);
            
            [x_vals,c_vals] = HNC_solve(params.Gamma,params.kappa);

            [coeff,scaling] = HNC_fit(x_vals,c_vals);
            %}
            
            weightsx = GaussLegendre_wgts*data.deltav1;
            weightsy = 1;
            weightsz = 1;
            if data.space_dims >= 2
                weightsy = GaussLegendre_wgts*data.deltav2;
                if data.space_dims >= 3
                    weightsz = GaussLegendre_wgts*data.deltav3;
                end
            end    
            

            
            for iq1 = 1:length(quadlocsx)
            for iq2 = 1:length(quadlocsy)
            for iq3 = 1:length(quadlocsz)
                
                summation = 0;
                for iqGH1 = 1:length(GaussHermite_locs)
                for iqGH2 = 1:length(GaussHermite_locs)
                for iqGH3 = 1:length(GaussHermite_locs)
                    quadx = quadlocsx(iq1) + GaussHermite_locs(iqGH1)/scaling;
                    quady = quadlocsy(iq2) + GaussHermite_locs(iqGH2)/scaling;
                    quadz = quadlocsz(iq3) + GaussHermite_locs(iqGH3)/scaling;
                    q = DGeval_3D(DGfluctuations,[quadx quady quadz],data); 
                    tilde_n_quad = q(1);
                    wgtx = GaussHermite_wgts(iqGH1);
                    wgty = GaussHermite_wgts(iqGH2);
                    wgtz = GaussHermite_wgts(iqGH3);                    
                    summation = summation + coeff/scaling*tilde_n_quad*wgtx*wgty*wgtz;                                    
                end
                end
                end
                
                wgtx = weightsx(iqGH1);
                wgty = weightsy(iqGH2);
                wgtz = weightsz(iqGH3);                    
                
                DFDN_cell = DFDN_cell + vectphi_Trans(1:data.Ls,1,iq1,iq2,iq3)*summation*wgtx*wgty*wgtz;
            end
            end
            end
            DFdn(:,iv1,iv2,iv3) = DFDN_cell;
            dpdn_excess_cell(1,iv1,iv2,iv3) = coeff;
            
            %done = (iv2+(iv1-1)*data.Nv2)/(data.Nv1*data.Nv2);
            %waitbar(done,h,'Computing Fcorr')                               
        end
    end
end
%delete(h)

