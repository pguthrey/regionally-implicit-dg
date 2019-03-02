function [ temp ] = DGeval_1D2V(DG1D2Vcoeffs,x,v1,v2,data)

%Periodic Conditions in x
domain = data.x_ub - data.x_lb;

xtilde = x - (ceil((x-data.x_lb)/domain)-1)*domain;
index = ceil((xtilde - data.x_lb)/data.deltax);
xii = 2*(xtilde-data.xcenters(index))/data.deltax;

v1cell = ceil((v1-data.v1_lb)/data.deltav1);
eta1 = (v1 - data.v1centers(v1cell))*2/data.deltav1;

v2cell = ceil((v2-data.v2_lb)/data.deltav2);
eta2 = (v2 - data.v2centers(v2cell))*2/data.deltav2;

%data.maxv2index = max([data.maxv2index v2cell])

zeta = zeros(1,data.theta1D2V);
for kay = 1:data.theta1D2V
    zeta(1,kay) = testfunction_zeta1D2V(xii,eta1,eta2,kay,data);
end



coeffs = DG1D2Vcoeffs(:,index,v1cell,v2cell);
temp = zeta*coeffs;