function [data] = finalize_error_save(DGsoln,data)
% written by Pierson Guthrey

solution = @(quadpoint) problem_solution(data.Tfinal,quadpoint,data);    

L1err_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L2err_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
Lierr_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L1exact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L2exact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
Liexact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);

data = initialize_gauss_lobatto_quadrature(data,data.M+2);
quadlocs = data.locs;
quadwgts = data.wgts1D;

plotexac= [];
plotappr = [];
plotx = [];

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

N_quad = length(quadlocs);
N_fullquad = N_quad.^data.space_dims;

N_quad_x = length(quadlocs);
N_quad_y = 1;
N_quad_z = 1;
if data.space_dims >= 2
    N_quad_y = N_quad;
    if data.space_dims >= 3
        N_quad_z = N_quad;
    end
end

filestring = ['problems/' data.problemname '/solutions/solution_' num2str(data.M) '_' data.correctorbasis '_' num2str(data.Nv1) '.mat'];
if exist(filestring, 'file')==2
    load(filestring,'exact','weightlist','vectphi_err');
else    
    %disp(' Installing error modules... warning: this will not work for systems yet ')
        
    %Make list of solution values at qudrature points
    exact = NaN(N_fullquad,Nv1,Nv2,Nv3);
    for iv1 = 1:data.Nv1
    for iv2 = 1:data.Nv2
    for iv3 = 1:data.Nv3
        v1bar = data.v1centers(iv1);
        v2bar = data.v1centers(iv2);
        v3bar = data.v1centers(iv3);

        quadlocsx = v1bar + data.deltav1/2*quadlocs;
        quadlocsy = NaN;
        quadlocsz = NaN;
        if data.space_dims >= 2
            quadlocsy = v2bar + data.deltav2/2*quadlocs;
            if data.space_dims >= 3
                quadlocsz = v3bar + data.deltav3/2*quadlocs;
            end
        end
        for iq1 = 1:N_quad_x
        for iq2 = 1:N_quad_y
        for iq3 = 1:N_quad_z
            quadpoint = [quadlocsx(iq1) quadlocsy(iq2) quadlocsz(iq3)];
            
            quadind = iq1 + (iq2-1)*N_quad_y ... 
                + (iq3-1)*N_quad_z*N_quad_y;
            
            exact(quadind,iv1,iv2,iv3) = solution(quadpoint); 
        end
        end
        end
    end
    end
    end    
    data.exact_save_compare = exact;
    
    weightlist = zeros(N_fullquad,1);
    weightsx = quadwgts*data.deltav1;
    weightsy = 1;
    weightsz = 1;
    if data.space_dims >= 2
        weightsy = quadwgts*data.deltav2;
        if data.space_dims >= 3
            weightsz = quadwgts*data.deltav3;
        end
    end        
    %Make list of quadrature weights
    for iq1 = 1:N_quad_x
    for iq2 = 1:N_quad_y
    for iq3 = 1:N_quad_z
        quadind = iq1 + (iq2-1)*N_quad_y ... 
            + (iq3-1)*N_quad_z*N_quad_y;

        weightlist(quadind,1) = ones(data.Neqns,1).*weightsx(iq1)*weightsy(iq2)*weightsz(iq3);   
    end
    end
    end
        
    data.weightlist_save_compare = weightlist;
    
    %Make conversion matrix from DG to sample pts
    testquad = zeros(data.M,N_quad);
    for iell = 1:data.M
    for iloc = 1:N_quad
        testquad(iell,iloc) = testfunction_Legendre(quadlocs(iloc),iell);    
    end
    end
    Zphi = 1;
    Yphi = 1;
    Xphi = 1;
    vectphi_err = zeros(N_fullquad,data.theta);
    for kay = 1:data.theta
        for eqn = 1:data.Neqns
        ells = data.sysBOs(eqn,kay);
        coeffs = data.sysBOscoeffs(eqn,kay);
            for v1quad = 1:N_quad_x
            for v2quad = 1:N_quad_y
            for v3quad = 1:N_quad_z
                if coeffs ~= 0 
                    Xphi = testquad(data.BOs(ells,1),v1quad);
                    if data.space_dims >= 2
                        Yphi = testquad(data.BOs(ells,2),v2quad);
                        if data.space_dims >= 3
                        Zphi = testquad(data.BOs(ells,3),v3quad);
                        end
                    end
                end

                quadind = v1quad + (v2quad-1)*N_quad_y ... 
                       + (v3quad-1)*N_quad_z*N_quad_y;
                   
                vectphi_err(quadind,kay) =  Xphi.*Yphi.*Zphi.*coeffs;
                data.vectphi_err_save_compare(quadind,kay) = vectphi_err(quadind,kay);
                
            end
            end
            end
        end
    end        
    save(filestring,'exact','weightlist','vectphi_err')
end


for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
for iv3 = 1:data.Nv3
    
        %{
    v1bar = data.v1centers(iv1);
    v2bar = data.v1centers(iv2);
    v3bar = data.v1centers(iv3);
    

    quadlocsx = v1bar + data.deltav1/2*quadlocs;
    quadlocsy = NaN;
    quadlocsz = NaN;
    if data.space_dims >= 2
        quadlocsy = v2bar + data.deltav2/2*quadlocs;
        if data.space_dims >= 3
            quadlocsz = v3bar + data.deltav3/2*quadlocs;
        end
    end
    weightsx = quadwgts*data.deltav1;
    weightsy = 1;
    weightsz = 1;
    if data.space_dims >= 2
        weightsy = quadwgts*data.deltav2;
        if data.space_dims >= 3
            weightsz = quadwgts*data.deltav3;
        end
    end    
    weight = NaN(data.Neqns,N_quad_x,N_quad_y,N_quad_z);
    approx = NaN(data.Neqns,N_quad_x,N_quad_y,N_quad_z);
    for iq1 = 1:N_quad_x
    for iq2 = 1:N_quad_y
    for iq3 = 1:N_quad_z
        quadpoint = [quadlocsx(iq1) quadlocsy(iq2) quadlocsz(iq3)];
        weight(:,iq1,iq2,iq3) = ones(data.Neqns,1).*weightsx(iq1)*weightsy(iq2)*weightsz(iq3);
        
        approx(:,iq1,iq2,iq3) = vectphi_err(:,:,iq1,iq2,iq3)*DGsoln(:,iv1,iv2,iv3);
        %approx(iq1,iq2,iq3) =  DGeval(DGsoln,quadpoint,data);        
        %truex_ind = iq1 + N_quad_x*(iv1-1);
        %truey_ind = iq2 + N_quad_y*(iv2-1);
        %truez_ind = iq3 + N_quad_z*(iv3-1);
        
        %{
        plotexac(truex_ind,truey_ind,truez_ind) = exact(:,iq1,iq2,iq3);
        plotappr(truex_ind,truey_ind,truez_ind) = approx(:,iq1,iq2,iq3);
        plotx(truex_ind,truey_ind,truez_ind) = quadlocsx(iq1);
        ploty(truex_ind,truey_ind,truez_ind) = quadlocsy(iq2);
        plotz(truex_ind,truey_ind,truez_ind) = quadlocsz(iq3);
        %}

    end
    end
    end
    %}
    
    approx_cell = vectphi_err*DGsoln(:,iv1,iv2,iv3);    
    error = abs(exact(:,iv1,iv2,iv3)-approx_cell);
    exactnorm = abs(exact(:,iv1,iv2,iv3));   
    for eqn = 1:data.Neqns
        
        L1err_norm = sum(error(:,1).*weightlist);
        L2err_norm = sum(error(:,1).^2.*weightlist);
        Lierr_norm = max(error(:,1));
        
        L1exact_norm = sum(exactnorm(:,1).*weightlist);
        L2exact_norm = sum(exactnorm(:,1).^2.*weightlist);
        Liexact_norm = max(exactnorm(:,1));
        
        L1err_norm_cell(eqn,iv1,iv2,iv3) = L1err_norm;
        L2err_norm_cell(eqn,iv1,iv2,iv3) = L2err_norm;
        Lierr_norm_cell(eqn,iv1,iv2,iv3) = Lierr_norm;
        
        L1exact_norm_cell(eqn,iv1,iv2,iv3) = L1exact_norm;
        L2exact_norm_cell(eqn,iv1,iv2,iv3) = L2exact_norm;
        Liexact_norm_cell(eqn,iv1,iv2,iv3) = Liexact_norm;
    end
end
end
end
data.L1err_norm_save_compare = L1err_norm;

for eqn = 1:data.Neqns
    L1err_norm_total(eqn) = sum(sum(sum(L1err_norm_cell(eqn,:,:,:))));
    L2err_norm_total(eqn) = sqrt(sum(sum(sum(L2err_norm_cell(eqn,:,:,:)))));
    Lierr_norm_total(eqn) = max(max(max(Lierr_norm_cell(eqn,:,:,:))));
    L1exact_norm_total(eqn) = sum(sum(sum(L1exact_norm_cell(eqn,:,:,:))));
    L2exact_norm_total(eqn) = sqrt(sum(sum(sum(L2exact_norm_cell(eqn,:,:,:)))));
    Liexact_norm_total(eqn) = max(max(max(Liexact_norm_cell(eqn,:,:,:))));
end

Rel_L1err = L1err_norm_total./L1exact_norm_total;
Rel_L2err = L2err_norm_total./L2exact_norm_total;
Rel_Lierr = Lierr_norm_total./Liexact_norm_total;

data.Rel_L1err = Rel_L1err;
data.Rel_L2err = Rel_L2err;
data.Rel_Lierr = Rel_Lierr;

data.Rel_L1err_filtered = NaN;
data.Rel_L2err_filtered = NaN;
data.Rel_Lierr_filtered = NaN;

