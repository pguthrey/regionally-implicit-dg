function [ zeta ] = testfunction_zeta1D(xii,ell,data)

zeta = 0;
if ell ~= 0
    ell_xii = data.BO1D(ell,1);
    phi_xii = testfunction_Legendre(xii,ell_xii);
    zeta = phi_xii;
end

end

