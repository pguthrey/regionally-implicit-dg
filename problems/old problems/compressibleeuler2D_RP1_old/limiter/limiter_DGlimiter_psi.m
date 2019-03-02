function [ DGcoeffs_limited ] = limiter_DGlimiter_psi(DGcoeffs,data)
%2D Limiter

DGcoeffs_limited = NaN(data.thetaT,data.Nv1,data.Nv2);

%Step 0
%DG_ghosts = problem_boundaryconditions_psi(DGcoeffs,data);
DG_ghosts = DGcoeffs;
%Step 1
q_max = zeros(data.Neqns,data.Nv1,data.Nv2);
q_min = zeros(data.Neqns,data.Nv1,data.Nv2);
q_ave = zeros(data.Neqns,data.Nv1,data.Nv2);
qcell = zeros(data.Neqns,data.P+2,data.P+2);
for iv1  = 1:(data.Nv1)
for iv2  = 1:(data.Nv2)
    q_min_cell = inf;
    q_max_cell = -inf;
    for tquad = 1:(data.P+2)
    for v1quad = 1:(data.P+2)
    for v2quad = 1:(data.P+2) 
        qcell_quad = data.vectpsi_limiter(:,:,tquad,v1quad,v2quad)*DG_ghosts(:,iv1,iv2);
        q_min_cell = min(q_min_cell,qcell_quad);
        q_max_cell = max(q_max_cell,qcell_quad);
    end
    end
    end
    q_ave(:,iv1,iv2) = data.vectpsi_average*DG_ghosts(:,iv1,iv2);
    q_min(:,iv1,iv2) = q_min_cell;
    q_max(:,iv1,iv2) = q_max_cell;    
end
end

%{
close all
A(:,:) = q_min;
surf(A)
axis([1 40 1 40 0 1])
figure
A(:,:) = q_ave;
surf(A)
axis([1 40 1 40 0 1])
figure
A(:,:) = q_max;
surf(A)
axis([1 40 1 40 0 1])
keyboard
%}

%Step 2
deltah = min([data.deltav1 data.deltav2]);
alpha = limiter_limiteralphafunction_psi(deltah)
data.M_max = zeros(data.Neqns,data.Nv1,data.Nv2);
data.M_min = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
for iv2  = 1:(data.Nv2)
    Qneighborsmax(:,1) = q_max(:,max([iv1-1 1]),max([iv2-1 1]));
    Qneighborsmax(:,2) = q_max(:,iv1-0,max([iv2-1 1]));
    Qneighborsmax(:,3) = q_max(:,min([iv1+1 data.Nv1]),max([iv2-1 1]));
    Qneighborsmax(:,4) = q_max(:,max([iv1-1 1]),iv2-0);
    Qneighborsmax(:,5) = q_ave(:,iv1,iv2)+alpha;
    Qneighborsmax(:,6) = q_max(:,min([iv1+1 data.Nv1]),iv2-0);
    Qneighborsmax(:,7) = q_max(:,max([iv1-1 1]),min([iv2+1 data.Nv2]));
    Qneighborsmax(:,8) = q_max(:,iv1-0,min([iv2+1 data.Nv2]));
    Qneighborsmax(:,9) = q_max(:,min([iv1+1 data.Nv1]),min([iv2+1 data.Nv2]));
    
    Qneighborsmin(:,1) = q_min(:,max([iv1-1 1]),max([iv2-1 1]));
    Qneighborsmin(:,2) = q_min(:,iv1-0,max([iv2-1 1]));
    Qneighborsmin(:,3) = q_min(:,min([iv1+1 data.Nv1]),max([iv2-1 1]));
    Qneighborsmin(:,4) = q_min(:,max([iv1-1 1]),iv2-0);
    Qneighborsmin(:,5) = q_ave(:,iv1,iv2)-alpha;
    Qneighborsmin(:,6) = q_min(:,min([iv1+1 data.Nv1]),iv2-0);
    Qneighborsmin(:,7) = q_min(:,max([iv1-1 1]),min([iv2+1 data.Nv2]));
    Qneighborsmin(:,8) = q_min(:,iv1-0,min([iv2+1 data.Nv2]));
    Qneighborsmin(:,9) = q_min(:,min([iv1+1 data.Nv1]),min([iv2+1 data.Nv2]));
    M_max(:,iv1,iv2) = max(Qneighborsmax,[],2);
    M_min(:,iv1,iv2) = min(Qneighborsmin,[],2);
end
end

%Step 3,4 and 5
cutoff_function = @(y) min(y./1.1,1);

limtol = 1e-8;
maxprecut = abs(M_max-q_ave)./(abs(q_max-q_ave)+limtol);
minprecut = abs(M_min-q_ave)./(abs(q_min-q_ave)+limtol);
%maxprecut = abs(1-q_ave)./(abs(q_max-q_ave)+limtol);
%minprecut = abs(q_ave)./(abs(q_min-q_ave)+limtol);

posprecut = abs(q_ave)./(abs(q_ave-q_min)+limtol);
inner_prod_average = data.vectpsi_average'*data.vectpsi_average;

%limiterpos = NaN(data.thetaT,data.Nv1,data.Nv2);

%Rossmanith/Seal Rescale
%%{
rescale = NaN(data.Nv1,data.Nv2);
for iv1  = 1:data.Nv1
for iv2  = 1:data.Nv2
    limitermax = cutoff_function(maxprecut(:,iv1,iv2));
    limitermin = cutoff_function(minprecut(:,iv1,iv2));
    limiterpos = posprecut(:,iv1,iv2);
    rescale(iv1,iv2) = max(min([ones(data.Neqns,1) limitermax limitermin limiterpos],[],2));
    %rescale(iv1,iv2) = max(min([ones(data.Neqns,1) limitermax limitermin],[],2));

    coeffs = DG_ghosts(:,iv1,iv2);

    qaverage = inner_prod_average*coeffs;
    
    DGtilde = (1-rescale(iv1,iv2)).*qaverage + rescale(iv1,iv2).*coeffs; 
    
    
    DGcoeffs_limited(:,iv1,iv2) = DGtilde;
end
end
%}
%Guthrey rescale
%{
N = 2^data.space_dims;
highermoments = (eye(data.thetaT)-inner_prod_average);
rescale = NaN(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:data.Nv1
for iv2  = 1:data.Nv2
    limitermax = cutoff_function(maxprecut(:,iv1,iv2));
    limitermin = cutoff_function(minprecut(:,iv1,iv2));
    limiterpos = posprecut(:,iv1,iv2);
    %escale(iv1,iv2) = max(min([ones(data.Neqns,1) limitermax limitermin limiterpos],[],2));
    vectbeta = min([ones(data.Neqns,1) limitermax limitermin],[],2);
    rescale(:,iv1,iv2) = vectbeta;
    B = diag(1-vectbeta);
    G = zeros(data.thetaT);
    for k = 1:data.Pd
        v1quad = data.Dlist(k,1); 
        if data.space_dims >= 2
            v2quad = data.Dlist(k,2);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
        else
            v3quad = 1;
        end
	weight = data.Dquadwgts(k);  
    psi = data.vectpsi(:,:,v1quad,v2quad,v3quad);
	G = G + weight*psi'*B*psi;        
    end	
    coeffs = DG_ghosts(:,iv1,iv2);
	DGtilde = coeffs+G*highermoments/N*coeffs;
    DGcoeffs_limited(:,iv1,iv2) = DGtilde;
end
end
%}
%{
close all
for k = 1:data.Neqns
    figure(k)
    R(:,:) = rescale(k,:,:);
    surf(R)
end
%}
%max(max(max(abs(DGcoeffs_limited(1,:,:)-DG_ghosts(1,:,:)))))
%keyboard



    