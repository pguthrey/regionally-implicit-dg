function [P_hartree_new] = compute_hartree(DGsolution,P_hartree_past,data)

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

deltav1 = data.deltav1;
deltav2 = data.deltav2;
deltav3 = data.deltav3;

N_total = Nv1*Nv2*Nv3;

%periodic = @(i,N) mod(i-1,N) + 1;
n_naut = zeros(N_total,1);
energyratio = zeros(N_total,1);
Upast = zeros(N_total,1);
for iv1 = 1:Nv1
    for iv2 = 1:Nv2
        for iv3 = 1:Nv3
            index = iv1 + (iv2-1)*Nv1 + (iv3-1)*Nv1*Nv2;
            Upast(index,1) = P_hartree_past(1,iv1,iv2,iv3);
            n_naut(index,1) = DGsolution(1,iv1,iv2,iv3);
            params = cell_average_to_params(n_naut(index,1),data);
            energyratio(index,1) = 1;%params.Eratio;
        end
    end
end


switch data.space_dims
    case 1
        LHS = @(U) (circshift(U,-1) - 2*U + circshift(U,+1))/deltav1^2;
    case 2
        LHS = @(U) (circshift(U,-1) - 2*U + circshift(U,+1))/deltav1^2 ...
            + (circshift(U,-Nv1) - 2*U + circshift(U,+Nv1))/deltav2^2;
    case 3
        LHS = @(U) (circshift(U,-1) - 2*U + circshift(U,+1))/deltav1^2 ... 
            + (circshift(U,-Nv1) - 2*U + circshift(U,+Nv1))/deltav2^2 ...
            + (circshift(U,-Nv1*Nv2) - 2*U + circshift(U,+Nv1*Nv2))/deltav3^2;   
end
                      
RHS = @(U) energyratio.*(n_naut - exp(-U));

options = optimoptions('fsolve','Display','off');

Unew = fsolve(@(U) LHS(U)-RHS(U), Upast , options);

P_hartree_new = P_hartree_past;
for iv1 = 1:Nv1
    for iv2 = 1:Nv2
        for iv3 = 1:Nv3
            index = iv1 + (iv2-1)*Nv1 + (iv3-1)*Nv1*Nv2;
            P_hartree_new(1,iv1,iv2,iv3) = Unew(index,1);
        end
    end
end








