function [ temp ] = DGeval_1D(coeffs,x,data)

xtilde = x - (ceil((x-data.x_lb)/(data.x_ub-data.x_lb))-1)*(data.x_ub-data.x_lb);
index = ceil((xtilde - data.x_lb)/data.deltax);
xii = 2*(xtilde-data.xcenters(index))/data.deltax;

temp = 0;
zeta = zeros(1,data.theta1D);
for kay = 1:data.theta1D
    zeta(1,kay) = testfunction_zeta1D(xii,kay,data);
end
temp = zeta*coeffs(:,index);
end

