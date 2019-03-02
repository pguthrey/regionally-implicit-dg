function [ data ] = check_conservation(DGcoeffs,data)
% written by Pierson Guthrey


Total = zeros(data.Neqns,1);

switch data.space_dims
    case 1
        dA = data.deltav1;
    case 2
        dA = data.deltav1*data.deltav2;
    case 3
        dA = data.deltav1*data.deltav2*data.deltav3;
end

for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
for iv3 = 1:data.Nv3    
    qcons = data.vectphi_average*DGcoeffs(:,iv1,iv2,iv3)*dA;
    Total = Total + qcons;
end
end
end

for eqn = 1:data.Neqns
    disp(['   Total ' data.varname{eqn} ' : ' num2str(Total(eqn))])
end

end

