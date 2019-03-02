function [Unew] = Electric_potential(DGsolution,Upast,data)

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

deltav1 = data.deltav1;
deltav2 = data.deltav2;
deltav3 = data.deltav3;

%periodic = @(i,N) mod(i-1,N) + 1;

for iv1 = 1:Nv1
    for iv2 = 1:Nv2
        for iv3 = 1:Nv3
            index = iv1 + (iv2-1)*Nv1 + (iv3-1)*Nv1*Nv2;
            n_naut(index,1) = DGsolution(1,:,:,:);
        end
    end
end



energyratio = data.energyratio;

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
                      
RHS = @(U) n_naut - exp(-U);

Unew = fsolve(@(U) LHS(U)-RHS(U), Upast );








