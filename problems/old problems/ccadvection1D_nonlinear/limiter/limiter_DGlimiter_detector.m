function [ BETA ] = limiter_DGlimiter_detector(DGcoeffs,DGpast,data)
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
        
        BETA(iv1) = max( BETA(iv1) , max(q_max(:,iv1) > M_max(:,iv1)+delta(iv1))); 
        BETA(iv1) = max( BETA(iv1) , max(q_cell < M_min(:,iv1)-delta(iv1)));     
    end
    end
end
