function [data] = problem_get_VMdata(data)

data.FinalTime = 45;
data.Nx = 120;
data.Neqns = 1;
data.nu = 1;
data.Nplotvars = 4;

data.vmboundsplot(1,:) = [-.01 .1];
data.vmboundsplot(2,:) = [-1.1 1.1];
data.vmboundsplot(3,:) = [-2 2];
data.vmboundsplot(4,:) = [-2 2];

%1D basis setup and testfunction precomputing
D = 1;
J=(1:data.M^(D+1))';
LM=zeros(data.M^(D+1),D+1); %cartesian product of bases
for j=1:D+1
    LM(:,D+1-j+1)=mod(ceil(J./data.M.^(j-1)),data.M);
end
LM(LM==0)=LM(LM==0)+data.M;
LO= LM-1; %Actual polynomial orders
data.BO= LM(sum(LO,2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
data.L1D = sum(data.BO(:,1)==1);
data.BO1D=data.BO(1:data.L1D,2:end);
data.theta1D = data.L1D*data.Neqns;

%1D basis setup and testfunction precomputing
D = 2;
J=(1:data.M^(D+1))';
LM=zeros(data.M^(D+1),D+1); %cartesian product of bases
for j=1:D+1
    LM(:,D+1-j+1)=mod(ceil(J./data.M.^(j-1)),data.M);
end
LM(LM==0)=LM(LM==0)+data.M;
LO= LM-1; %Actual polynomial orders
data.BO= LM(sum(LO,2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
data.L1D1V = sum(data.BO(:,1)==1);
data.BO1D1V =data.BO(1:data.L1D1V,2:end);
data.theta1D1V = data.L1D1V*data.Neqns;
%1V basis setup and testfunction precomputing


[data] = initialize_Basis(data,1);

[data] = initialize_gauss_quadrature(data,data.M);

%Derived quanities
data.Neqns = 1;
xendpts = linspace(data.x_lb,data.x_ub,data.Nx+1);
data.xcenters = (xendpts(1:end-1)+xendpts(2:end))/2;
data.deltax = (data.x_ub-data.x_lb)/data.Nx;
v1endpts = linspace(data.v1_lb,data.v1_ub,data.Nv1+1);
data.v1centers = (v1endpts(1:end-1)+v1endpts(2:end))/2;
data.deltav1 = (data.v1_ub-data.v1_lb)/data.Nv1;
deltat = data.deltav1/data.nu/data.Fspeedmax*data.cfl;
data.nuv1 = data.nu*deltat/data.deltav1;

data = initialize_precompute_testfunctions(data);
data = precompute_VM_testfunctions(data);
data = initialize_formIntegrals(data);
data = initialize_method(data);

P = data.P;
for xquad = 1:P
    for v1quad = 1:P
        data.wgts2D(xquad,v1quad) = data.wgts1D(xquad)*data.wgts1D(v1quad);
    end
end


end

