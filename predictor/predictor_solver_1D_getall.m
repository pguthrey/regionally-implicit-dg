function [soln,data,Jacobian,residual] = predictor_solver_1D_getall(this,data)

thetaT = data.thetaT;
cells_per_region = data.cells_per_region;
region_per_dimension = data.region_per_dimension;
rx_param = data.rx_param;
ry_param = data.ry_param;
rz_param = data.rz_param;
space_dims = data.space_dims;
main_cell = data.main_cell;
                                        
index = @(gamma) (gamma-1)*thetaT+(1:thetaT);
N = cells_per_region*thetaT;
%M = (cells_per_region*3-2)*thetaT^2;
Jacobian = zeros(N,N);
residual = zeros(N,1);
%Jacobian = spalloc(N,N,M);

cellnum = 1;
regions = region_per_dimension;
rz = rz_param;
ry = ry_param;
rx = rx_param;

for region_k = -rz:rz
for region_j = -ry:ry
for region_i = -rx:rx
    row = index(cellnum);
    col = index(cellnum);
%    cellj = iv2+region_j;
%    cellk = iv3+region_k;
    celli = region_i+1+rx;
    cellj = region_j+1+ry;
    cellk = region_k+1+rz;
    residual(row,1) = this.residual_cell(:,celli,cellj,cellk);
    Jacobian(row,col) = this.Jac_cell(:,:,celli,cellj,cellk);
    %West Jacobian flux
    if region_i > -rx
        residual(row,1) = residual(row,1) + this.west_flux(:,celli,cellj,cellk);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_west_cell(:,:,celli,cellj,cellk);
        col = index(cellnum-1);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_west_other(:,:,celli,cellj,cellk);
    else
        residual(row,1) = residual(row,1) + this.west_trunc(:,celli,cellj,cellk);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_west_trunc(:,:,celli,cellj,cellk);
    end
    col = index(cellnum);
    %East Jacobian flux
    if region_i < rx
        residual(row,1) = residual(row,1) + this.east_flux(:,celli,cellj,cellk);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_east_cell(:,:,celli,cellj,cellk);
        col = index(cellnum+1);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_east_other(:,:,celli,cellj,cellk);
    else
        residual(row,1) = residual(row,1) + this.east_trunc(:,celli,cellj,cellk);
        Jacobian(row,col) = Jacobian(row,col) + this.Jac_east_trunc(:,:,celli,cellj,cellk);
    end            
    col = index(cellnum);

    if space_dims >= 2
        %sout Jacobian flux
        if region_j > -ry
            residual(row,1) = residual(row,1) + this.sout_flux(:,celli,cellj,cellk);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_sout_cell(:,:,celli,cellj,cellk);
            col = index(cellnum-regions);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_sout_other(:,:,celli,cellj,cellk);
        else
            residual(row,1) = residual(row,1) + this.sout_trunc(:,celli,cellj,cellk);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_sout_trunc(:,:,celli,cellj,cellk);
        end
        col = index(cellnum);
        %nort Jacobian flux
        if region_j < ry
            residual(row,1) = residual(row,1) + this.nort_flux(:,celli,cellj,cellk);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_nort_cell(:,:,celli,cellj,cellk);
            col = index(cellnum+regions);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_nort_other(:,:,celli,cellj,cellk);
        else
            residual(row,1) = residual(row,1) + this.nort_trunc(:,celli,cellj,cellk);
            Jacobian(row,col) = Jacobian(row,col) + this.Jac_nort_trunc(:,:,celli,cellj,cellk);
        end            
        if space_dims >= 3        
            col = index(cellnum);
            %down Jacobian flux
            if region_k > -rz
                residual(row,1) = residual(row,1) + this.flux_down(:,celli,cellj,cellk);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_down_cell(:,:,celli,cellj,cellk);
                col = index(cellnum-regions^2);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_down_other(:,:,celli,cellj,cellk);
            else
                residual(row,1) = residual(row,1) + this.trunc_down(:,celli,cellj,cellk);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_trunc_down(:,:,celli,cellj,cellk);
            end
            col = index(cellnum);
            %uppr Jacobian flux
            if region_k < rz
                residual(row,1) = residual(row,1) + this.flux_uppr(:,celli,cellj,cellk);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_uppr_cell(:,:,celli,cellj,cellk);
                col = index(cellnum+regions^2);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_uppr_other(:,:,celli,cellj,cellk);
            else
                residual(row,1) = residual(row,1) + this.trunc_uppr(:,celli,cellj,cellk);
                Jacobian(row,col) = Jacobian(row,col) + this.Jac_trunc_uppr(:,:,celli,cellj,cellk);
            end  
        end
    end
    %iterate cell number
    cellnum = cellnum + 1;
end
end
end


if any(any(any(isnan(residual))))
    keyboard
end

soln = Jacobian\residual;


end

