function [ DGcoeffs_limited ] = limiter_DGlimiter(DGcoeffs,DGpast,data)
%2D Limiter

newlocs = [-1;data.locs;1];

DGcoeffs_limited = NaN(data.theta,data.Nv1,data.Nv2);

%Step 0
%DG_ghosts = problem_boundaryconditions_phi(DGcoeffs,data);
%Step 1
q_max = zeros(data.Neqns,data.Nv1,data.Nv2);
q_min = zeros(data.Neqns,data.Nv1,data.Nv2);
q_ave = zeros(data.Neqns,data.Nv1,data.Nv2);
qcell = zeros(data.Neqns,data.P+2);
for iv1  = 1:(data.Nv1)
    q_min_cell = inf;
    q_max_cell = -inf;
    for v1quad = 1:(data.P+2)
        qcell_quad = data.vectphi_limiter(:,:,v1quad)*DGcoeffs(:,iv1);
        q_min_cell = min(q_min_cell,qcell_quad);
        q_max_cell = max(q_max_cell,qcell_quad);
    end
    q_ave(:,iv1) = data.vectphi_average*DGcoeffs(:,iv1);
    q_min(:,iv1) = q_min_cell;
    q_max(:,iv1) = q_max_cell; 
end

q_max_past = zeros(data.Neqns,data.Nv1,data.Nv2);
q_min_past = zeros(data.Neqns,data.Nv1,data.Nv2);
q_ave_past = zeros(data.Neqns,data.Nv1,data.Nv2);
qcell_past = zeros(data.Neqns,data.P+2);
for iv1  = 1:(data.Nv1)
    q_min_cell_past = inf;
    q_max_cell_past = -inf;
    for v1quad = 1:(data.P+2)
        qcell_quad = data.vectphi_limiter(:,:,v1quad)*DGcoeffs(:,iv1);
        q_min_cell_past = min(q_min_cell_past,qcell_quad);
        q_max_cell_past = max(q_max_cell_past,qcell_quad);
    end
    q_ave_past(:,iv1) = data.vectphi_average*DGpast(:,iv1);
    q_min_past(:,iv1) = q_min_cell_past;
    q_max_past(:,iv1) = q_max_cell_past;    
end


%Step 2
deltah = data.deltav1;
alpha = limiter_limiteralphafunction_phi(deltah);
M_max = zeros(data.Neqns,data.Nv1,data.Nv2);
M_min = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
    Qneighborsmax(:,1) = q_max(:,max([iv1-1 1]));
    Qneighborsmax(:,2) = q_ave(:,iv1)+alpha;
    Qneighborsmax(:,3) = q_max(:,min([iv1+1 data.Nv1]));
    %%{
    Qneighborsmax(:,4) = q_max_past(:,max([iv1-1 1]));
    Qneighborsmax(:,5) = q_max_past(:,iv1);
    Qneighborsmax(:,6) = q_max_past(:,min([iv1+1 data.Nv1]));
    %}
    Qneighborsmin(:,1) = q_min(:,max([iv1-1 1]));
    Qneighborsmin(:,2) = q_ave(:,iv1)-alpha;
    Qneighborsmin(:,3) = q_min(:,min([iv1+1 data.Nv1]));
    %%{
    Qneighborsmin(:,4) = q_min_past(:,max([iv1-1 1]));
    Qneighborsmin(:,5) = q_min_past(:,iv1);
    Qneighborsmin(:,6) = q_min_past(:,min([iv1+1 data.Nv1]));    
    %}
    M_max(:,iv1) = max(Qneighborsmax,[],2);
    M_min(:,iv1) = min(Qneighborsmin,[],2);
end

%Step 3,4 and 5
cutoff_function = @(y) min(y./1.1,1);

limtol = min([1e-8 alpha/2]);
maxprecut = (abs(M_max-q_ave)+2*limtol)./(abs(q_max-q_ave)+limtol);
minprecut = (abs(M_min-q_ave)+2*limtol)./(abs(q_min-q_ave)+limtol);
%maxprecut = abs(1-q_ave)./(abs(q_max-q_ave)+limtol);
%minprecut = abs(q_ave)./(abs(q_min-q_ave)+limtol);

%limiterpos = NaN(data.theta,data.Nv1,data.Nv2);

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

%Rossmanith/Seal Rescale
%%{
rescale = NaN(data.Nv1,data.Nv2);
for iv1  = 1:data.Nv1
    limitermax(:,iv1) = cutoff_function(maxprecut(:,iv1));
    limitermin(:,iv1) = cutoff_function(minprecut(:,iv1));
    rescale(iv1) = max(min([ones(data.Neqns,1) limitermax(:,iv1) limitermin(:,iv1)],[],2));
    %rescale(iv1) = max(min([ones(data.Neqns,1) limitermax limitermin],[],2));

    %rescale(iv1) = rescale(iv1)>= 1;
    coeffs = DGcoeffs(:,iv1);

    qaverage = data.inner_prod_average_phi*coeffs;
    
    DGtilde = (1-rescale(iv1)).*qaverage + rescale(iv1).*coeffs; 
    
    DGcoeffs_limited(:,iv1) = DGtilde;
end

%{
figure(77+3)
plot(rescale)
keyboard
disp('eeeeeeeeeeeeerrrrrrrrrrr')
%}


    