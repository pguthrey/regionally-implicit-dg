function [ temp ] = DGeval_1D1V(DG1D2Vcoeffs,x,v1,data)

%Periodic Conditions in x
domain = data.x_ub - data.x_lb;

xtilde = x - (ceil((x-data.x_lb)/domain)-1)*domain;
index = ceil((xtilde - data.x_lb)/data.deltax);
xii = 2*(xtilde-data.xcenters(index))/data.deltax;

v1cell = ceil((v1-data.v1_lb)/data.deltav1);
eta1 = (v1 - data.v1centers(v1cell))*2/data.deltav1;

zeta = zeros(1,data.theta1D1V);
for kay = 1:data.theta1D1V
    zeta(1,kay) = testfunction_zeta1D1V(xii,eta1,kay,data);
end

coeffs = DG1D2Vcoeffs(:,index,v1cell);
temp = zeta*coeffs;