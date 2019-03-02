function [data] = problem_get_VMdata(data)

data.FinalTime = 45;
data.Nx = 101;
data.cfl = .5;
data.Neqns = 1;
data.nu = 1;
data.Nplotvars = 4;

data.vmboundsplot(1,:) = [-.05 .1];
data.vmboundsplot(2,:) = [-10 10];
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
D = 3;
J=(1:data.M^(D+1))';
LM=zeros(data.M^(D+1),D+1); %cartesian product of bases
for j=1:D+1
    LM(:,D+1-j+1)=mod(ceil(J./data.M.^(j-1)),data.M);
end
LM(LM==0)=LM(LM==0)+data.M;
LO= LM-1; %Actual polynomial orders
data.BO= LM(sum(LO,2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
data.L1D2V = sum(data.BO(:,1)==1);
data.BO1D2V =data.BO(1:data.L1D2V,2:end);
data.theta1D2V = data.L1D2V*data.Neqns;
%2V basis setup and testfunction precomputing
[data] = initialize_Basis(data,2);

data.x_ub =  2*pi;
data.x_lb = -2*pi;

[data] = initialize_gauss_quadrature(data);

%Derived quanities
data.Neqns = 1;
xendpts = linspace(data.x_lb,data.x_ub,data.Nx+1);
data.xcenters = (xendpts(1:end-1)+xendpts(2:end))/2;
data.deltax = (data.x_ub-data.x_lb)/data.Nx;
v1endpts = linspace(data.v1_lb,data.v1_ub,data.Nv1+1);
data.v1centers = (v1endpts(1:end-1)+v1endpts(2:end))/2;
data.deltav1 = (data.v1_ub-data.v1_lb)/data.Nv1;
v2endpts = linspace(data.v2_lb,data.v2_ub,data.Nv2+1);
data.v2centers = (v2endpts(1:end-1)+v2endpts(2:end))/2;
data.deltav2 = (data.v2_ub-data.v2_lb)/data.Nv2;
deltat = min([data.deltav1/data.nu/data.Fspeedmax*data.cfl data.deltav2/1/data.Gspeedmax*data.cfl]);
data.nuv1 = data.nu*deltat/data.deltav1;
data.nuv2 = 1*deltat/data.deltav2;

data = initialize_precompute_testfunctions(data);
data = precompute_VM_testfunctions(data);
data = initialize_formIntegrals(data);

P = data.P;
for xquad = 1:P
    for v1quad = 1:P
        for v2quad = 1:P
            data.wgts3D(xquad,v1quad,v2quad) = data.wgts1D(xquad)*data.wgts1D(v1quad)*data.wgts1D(v2quad);
        end
    end
end


end

