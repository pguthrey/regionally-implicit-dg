function [ zeta ] = testfunction_zeta1D1V(xii,eta1,ell,data)

zeta = 0;
if ell ~= 0
    ell_xii = data.BO1D1V(ell,1);
    ell_eta1 = data.BO1D1V(ell,2);
    phi_xii = testfunction_Legendre(xii,ell_xii);
    phi_eta1 = testfunction_Legendre(eta1,ell_eta1);
    zeta = phi_xii.*phi_eta1;
end

end

