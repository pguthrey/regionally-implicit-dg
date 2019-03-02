function [ temp ] = DGeval_2D(DGcoeffs,vpoint,data)
% written by Pierson Guthrey

v1 = vpoint(1);
v2 = vpoint(2);

%Periodic Conditions in x
domainv1 = data.v1_ub - data.v1_lb;
domainv2 = data.v2_ub - data.v2_lb;

v1tilde = v1 - (ceil((v1-data.v1_lb)/domainv1)-1)*domainv1;
indexv1 = ceil((v1tilde - data.v1_lb)/data.deltav1);
etav1 = 2*(v1tilde-data.v1centers(indexv1))/data.deltav1;

v2tilde = v2 - (ceil((v2-data.v2_lb)/domainv2)-1)*domainv2;
indexv2 = ceil((v2tilde - data.v2_lb)/data.deltav2);
etav2 = 2*(v2tilde-data.v2centers(indexv2))/data.deltav2;

zeta = zeros(1,data.theta);
quadpoint = [etav1,etav2];
for kay = 1:data.theta
    for eqn = 1:data.Neqns
        ell = data.sysBOs(eqn,kay);
        alpha = data.sysBOscoeffs(eqn,kay);
        if all([ell alpha]~=0)
            zeta(1,kay) = alpha*testfunction_phi(quadpoint,ell,data);
        end
    end
end

coeffs = DGcoeffs(:,indexv1,indexv2);
temp = zeta*coeffs;


