function [ DGcoeffs_limited ] = limiter_DGlimiter_psi(DGcoeffs,DGpast,data)
%2D Limiter

DGcoeffs_limited = NaN(data.thetaT,data.Nv1,data.Nv2);

%Step 0
%DG_ghosts = problem_boundaryconditions_psi(DGcoeffs,data);
%Step 1
q_max = zeros(data.Neqns,data.Nv1,data.Nv2);
q_min = zeros(data.Neqns,data.Nv1,data.Nv2);
q_ave = zeros(data.Neqns,data.Nv1,data.Nv2);
qcell = zeros(data.Neqns,data.P+2);
for iv1  = 1:(data.Nv1)
    q_min_cell = inf;
    q_max_cell = -inf;
    for tquad = 1:(data.P+2)
    for v1quad = 1:(data.P+2)
        qcell_quad = data.vectpsi_limiter(:,:,tquad,v1quad)*DGcoeffs(:,iv1);
        q_min_cell = min(q_min_cell,qcell_quad);
        q_max_cell = max(q_max_cell,qcell_quad);
    end
    end
    q_ave(:,iv1) = data.vectpsi_average*DGcoeffs(:,iv1);
    q_min(:,iv1) = q_min_cell;
    q_max(:,iv1) = q_max_cell;    
end

%Step 2
deltah = data.deltav1;
alpha = limiter_limiteralphafunction_psi(deltah);
M_max = zeros(data.Neqns,data.Nv1,data.Nv2);
M_min = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
    Qneighborsmax(:,1) = q_max(:,max([iv1-1 1]));
    Qneighborsmax(:,2) = q_ave(:,iv1)+alpha;
    Qneighborsmax(:,3) = q_max(:,min([iv1+1 data.Nv1]));
    Qneighborsmin(:,1) = q_min(:,max([iv1-1 1]));
    Qneighborsmin(:,2) = q_ave(:,iv1)-alpha;
    Qneighborsmin(:,3) = q_min(:,min([iv1+1 data.Nv1]));
    M_max(:,iv1) = max(Qneighborsmax,[],2);
    M_min(:,iv1) = min(Qneighborsmin,[],2);
end

%Step 3,4 and 5
cutoff_function = @(y) min(y./1.1,1);


limtol = 1e-8;
maxprecut = (abs(M_max-q_ave)+2*limtol)./(abs(q_max-q_ave)+limtol);
minprecut = (abs(M_min-q_ave)+2*limtol)./(abs(q_min-q_ave)+limtol);
%maxprecut = abs(1-q_ave)./(abs(q_max-q_ave)+limtol);
%minprecut = abs(q_ave)./(abs(q_min-q_ave)+limtol);

%{
figure(77)
plot(minprecut)
title('minprecut')
axis([1 data.Nv1 0 4])
figure(77+1)
plot(maxprecut)
title('maxprecut')
axis([1 data.Nv1 0 4])

figure(77+2)
temp = [maxprecut ; minprecut];
plot(min(temp))
title('precut')
axis([1 data.Nv1 0 4]) 
debug_plot_DGpsi(DGcoeffs,data)
%}

%Rossmanith/Seal Rescale
%%{
rescale = NaN(data.Nv1,data.Nv2);
for iv1  = 1:data.Nv1
    limitermax = cutoff_function(maxprecut(:,iv1));
    limitermin = cutoff_function(minprecut(:,iv1));
    %limiterpos = posprecut(:,iv1);
    rescale(iv1) = max(min([ones(data.Neqns,1) limitermax limitermin],[],2));

    %rescale(iv1) = rescale(iv1)>= 1;
    coeffs = DGcoeffs(:,iv1);

    qaverage = data.inner_prod_average_psi*coeffs;
    
    DGtilde = (1-rescale(iv1)).*qaverage + rescale(iv1).*coeffs; 
        
    DGcoeffs_limited(:,iv1) = DGtilde;
       
end

%debug_plot_DGpsi(DGcoeffs_limited,data)

DGcoeffs_limited = DGcoeffs;
%keyboard
