function [ DGcoeffs_limited ] = limiter_DGlimiter_psi(DGcoeffs,DGpast,data)


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


%Positivity Violation Detector
BETA = zeros(1,data.Nv1);
BETA = max([ BETA ; max( q_min < 0 , [] , 1) ] , [] , 1);

q_max_past = zeros(data.Neqns,data.Nv1,data.Nv2);
q_min_past = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
    q_min_cell_past = inf;
    q_max_cell_past = -inf;
    for v1quad = 1:(data.P+2)
        qcell_quad = data.vectphi_limiter(:,:,v1quad)*DGpast(:,iv1);
        q_min_cell_past = min(q_min_cell_past,qcell_quad);
        q_max_cell_past = max(q_max_cell_past,qcell_quad);
    end
    q_min_past(:,iv1) = q_min_cell_past;
    q_max_past(:,iv1) = q_max_cell_past;    
end

%{
%Step 2
M_max = zeros(data.Neqns,data.Nv1,data.Nv2);
M_min = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
    Qneighborsmax(:,1) = q_max(:,max([iv1-1 1]));
    Qneighborsmax(:,2) = q_max(:,iv1);
    Qneighborsmax(:,3) = q_max(:,min([iv1+1 data.Nv1]));
    %%{
    Qneighborsmax(:,4) = q_max_past(:,max([iv1-1 1]));
    Qneighborsmax(:,5) = q_max_past(:,iv1);
    Qneighborsmax(:,6) = q_max_past(:,min([iv1+1 data.Nv1]));
    Qneighborsmin(:,1) = q_min(:,max([iv1-1 1]));
    Qneighborsmin(:,2) = q_min(:,iv1);
    Qneighborsmin(:,3) = q_min(:,min([iv1+1 data.Nv1]));
    %%{
    Qneighborsmin(:,4) = q_min_past(:,max([iv1-1 1]));
    Qneighborsmin(:,5) = q_min_past(:,iv1);
    Qneighborsmin(:,6) = q_min_past(:,min([iv1+1 data.Nv1]));    
    M_max(:,iv1) = max(Qneighborsmax,[],2);
    M_min(:,iv1) = min(Qneighborsmin,[],2);
end
%}

%Step 2
M_max = zeros(data.Neqns,data.Nv1,data.Nv2);
M_min = zeros(data.Neqns,data.Nv1,data.Nv2);
for iv1  = 1:(data.Nv1)
    Qneighborsmax(:,1) = q_max_past(:,max([iv1-1 1]));
    Qneighborsmax(:,2) = q_max_past(:,iv1);
    Qneighborsmax(:,3) = q_max_past(:,min([iv1+1 data.Nv1]));
    
    Qneighborsmin(:,1) = q_min_past(:,max([iv1-1 1]));
    Qneighborsmin(:,2) = q_min_past(:,iv1);
    Qneighborsmin(:,3) = q_min_past(:,min([iv1+1 data.Nv1]));
    
    M_max(:,iv1) = max(Qneighborsmax,[],2);
    M_min(:,iv1) = min(Qneighborsmin,[],2);
end


delta = min(1e-4,(M_max-M_min)*1e-3 + 1e-8);

for iv1  = 1:(data.Nv1)
    for tquad = 1:(data.P+2)
    for v1quad = 1:(data.P+2)        
        q_cell = data.vectpsi_limiter(:,:,tquad,v1quad)*DGcoeffs(:,iv1);
        
        BETA_max(iv1) = max( BETA(iv1) , max(q_max(:,iv1) > M_max(:,iv1)+delta(iv1))); 
        BETA_min(iv1) = max( BETA(iv1) , max(q_cell < M_min(:,iv1)-delta(iv1)));     
    end
    end
end


limtol = 1e-8;

xmesh = 1:data.Nv1;
DGcoeffs_limited = DGcoeffs;
for iv1  = xmesh(BETA==1)
    Averages = data.V*DGcoeffs_limited(:,iv1)/data.limiter_dc^2;
    f_ave = sum(Averages)/(2*data.M+1)^(data.space_dims+1);
    f_max = max(Averages);
    f_min = min(Averages);
    
    maxprecut = (abs(M_max-f_ave)+2*limtol)./(abs(f_max-q_ave(:,iv1))+limtol);
    minprecut = (abs(M_min-f_ave)+2*limtol)./(abs(f_min-q_ave(:,iv1))+limtol);
    
    precut = min(min([maxprecut minprecut])/1.1, 1);
    
    Averages2 = Averages*precut+(1-precut)*f_ave;
    
    Psi_cell = data.V'*Averages2/2^(data.space_dims+1);
    DGcoeffs_limited(:,iv1) = Psi_cell;
end

%sDGcoeffs_limited = DGcoeffs;
%{
    debug_plot_DGpsi(DGcoeffs_limited,data)    
    debug_plot_DGpsi(DGcoeffs,data)
%}
%{


DGcoeffs_limited = DGcoeffs;
debug_plot_DGpsi(DGcoeffs_limited,data)
disp('pre limiting')
pause

xmesh = 1:data.Nv1;
DGcoeffs_limited = DGcoeffs;
check = 1;
iters = 0;
[ BETA ] = limiter_DGlimiter_detector(DGcoeffs_limited,DGpast,data);
while check

    iters = iters + 1
    
    for iv1  = xmesh(BETA==1)
        Averages = data.V*DGcoeffs_limited(:,iv1)
        
        
        DGcoeffs_limited(:,iv1) = data.limiter_filter*;
    end
    
    [ BETA ] = limiter_DGlimiter_detector(DGcoeffs_limited,DGpast,data);




    check = any(BETA == 1) && (iters < 1000);
end


keyboard
%}