function [ temp ] = DGeval_1D(DGcoeffs,vpoint,data)
% written by Pierson Guthrey

v1 = vpoint(1);

%Periodic Conditions in x
domainv1 = data.v1_ub - data.v1_lb;
v1tilde = -mod(data.v1_ub - v1,domainv1)+data.v1_ub;

indexv1 = ceil((v1tilde - data.v1_lb)/data.deltav1);
etav1 = 2*(v1tilde-data.v1centers(indexv1))/data.deltav1;

zeta = zeros(1,data.theta);
quadpoint = etav1;
for kay = 1:data.theta
    for eqn = 1:data.Neqns
        ell = data.sysBOs(eqn,kay);
        alpha = data.sysBOscoeffs(eqn,kay);
        if all([ell alpha]~=0)
            zeta(1,kay) = alpha*testfunction_phi(quadpoint,ell,data);
        end
    end
end

coeffs = DGcoeffs(:,indexv1);
temp = zeta*coeffs;


