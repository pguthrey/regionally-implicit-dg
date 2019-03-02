function [data] = finalize_error(DGsoln,data)
% written by Pierson Guthrey

solution = @(quadpoint) problem_solution(data.Tfinal,quadpoint,data);    

L1err_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L2err_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
Lierr_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L1exact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
L2exact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);
Liexact_norm_cell = NaN(data.Neqns,data.Nv1,data.Nv2,data.Nv3);

data = initialize_gauss_lobatto_quadrature(data,data.M+3);
quadlocs = data.locs;
quadwgts = data.wgts1D;

plotexac= [];
plotappr = [];
plotx = [];

testquad = zeros(1,data.M,length(quadlocs));
for iell = 1:data.M
for iloc = 1:length(quadlocs)
    testquad(1,iell,iloc) = testfunction_Legendre(quadlocs(iloc),iell);    
end
end

Zphi = 1;
Yphi = 1;
Xphi = 1;
vectphi_err = zeros(1,data.theta,length(quadlocs),length(quadlocs),length(quadlocs));
for kay = 1:data.theta
    for eqn = 1:data.Neqns
    ells = data.sysBOs(eqn,kay);
    coeffs = data.sysBOscoeffs(eqn,kay);
        for v1quad = 1:length(quadlocs)
        for v2quad = 1:length(quadlocs)
        for v3quad = 1:length(quadlocs)
            if coeffs ~= 0 
                Xphi = testquad(1,data.BOs(ells,1),v1quad);
                if data.space_dims >= 2
                    Yphi = testquad(1,data.BOs(ells,2),v2quad);
                    if data.space_dims >= 3
                    Zphi = testquad(1,data.BOs(ells,3),v3quad);
                    end
                end
            end
        vectphi_err(eqn,kay,v1quad,v2quad,v3quad) =  Xphi.*Yphi.*Zphi.*coeffs;
        
        quadind = v1quad + (v2quad-1)*length(quadlocs);
        data.vectphi_err_compare(quadind,kay) = Xphi.*Yphi.*Zphi.*coeffs;
        
        end
        end
        end
    end
end

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
    weightsx = quadwgts*data.deltav1;
    weightsy = 1;
    weightsz = 1;
    if data.space_dims >= 2
        weightsy = quadwgts*data.deltav2;
        if data.space_dims >= 3
            weightsz = quadwgts*data.deltav3;
        end
    end    

    weight = NaN(data.Neqns,length(quadlocsx),length(quadlocsy),length(quadlocsz));
    approx = NaN(data.Neqns,length(quadlocsx),length(quadlocsy),length(quadlocsz));
    exact = NaN(data.Neqns,length(quadlocsx),length(quadlocsy),length(quadlocsz));
    for iq1 = 1:length(quadlocsx)
    for iq2 = 1:length(quadlocsy)
    for iq3 = 1:length(quadlocsz)
        quadpoint = [quadlocsx(iq1) quadlocsy(iq2) quadlocsz(iq3)];
        weight(:,iq1,iq2,iq3) = ones(data.Neqns,1).*weightsx(iq1)*weightsy(iq2)*weightsz(iq3);
        approx(:,iq1,iq2,iq3) = vectphi_err(:,:,iq1,iq2,iq3)*DGsoln(:,iv1,iv2,iv3);
        %approx(iq1,iq2,iq3) =  DGeval(DGsoln,quadpoint,data);
        exact(:,iq1,iq2,iq3) = solution(quadpoint);  
        
        truex_ind = iq1 + length(quadlocsx)*(iv1-1);
        truey_ind = iq2 + length(quadlocsy)*(iv2-1);
        truez_ind = iq3 + length(quadlocsz)*(iv3-1);
        
        %{
        plotexac(truex_ind,truey_ind,truez_ind) = exact(:,iq1,iq2,iq3);
        plotappr(truex_ind,truey_ind,truez_ind) = approx(:,iq1,iq2,iq3);
        plotx(truex_ind,truey_ind,truez_ind) = quadlocsx(iq1);
        ploty(truex_ind,truey_ind,truez_ind) = quadlocsy(iq2);
        plotz(truex_ind,truey_ind,truez_ind) = quadlocsz(iq3);
        %}
        quadind = iq1 + (iq2-1)*length(quadlocs);
        data.weightlist_compare(quadind,1) = weight(:,iq1,iq2,iq3);

        data.exact_compare(quadind,iv1,iv2) = exact(:,iq1,iq2,iq3);

    end
    end
    end
    error = abs(exact-approx);
    exactnorm = abs(exact);
   
    for eqn = 1:data.Neqns
        L1err_norm = sum(sum(sum(error(eqn).*weight)));
        L2err_norm = sum(sum(sum(error(eqn).^2.*weight)));
        Lierr_norm = max(max(max(error(eqn))));
        L1exact_norm = sum(sum(sum(exactnorm(eqn).*weight)));
        L2exact_norm = sum(sum(sum(exactnorm(eqn).^2.*weight)));
        Liexact_norm = max(max(max(exactnorm(eqn))));
        
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
data.L1err_norm_compare = L1err_norm;


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


%{

N = length(quadlocs)*data.Nv2;

close all
subplot(1,2,1) 
mixed = [plotappr(1:floor(N/2),: ); plotexac( floor(N/2)+1 : end, : )];
surf(plotx,ploty,mixed) 
title('approximate exact ')
shading interp 
view(2) 
colorbar 

subplot(1,2,2) 
surf(plotx,ploty,abs(plotexac-plotappr)) 
title('diff')
shading interp 
view(2) 
colorbar 

keyboard
%}

data.Rel_L1err = Rel_L1err;
data.Rel_L2err = Rel_L2err;
data.Rel_Lierr = Rel_Lierr;

data.Rel_L1err_filtered = NaN;
data.Rel_L2err_filtered = NaN;
data.Rel_Lierr_filtered = NaN;

