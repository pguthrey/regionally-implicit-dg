function [ zeta ] = testfunction_zeta1D2V(xii,eta1,eta2,ell,data)

zeta = 0;
if ell ~= 0
    ell_xii = data.BO1D2V(ell,1);
    ell_eta1 = data.BO1D2V(ell,2);
    ell_eta2 = data.BO1D2V(ell,3);
    phi_xii = testfunction_Legendre(xii,ell_xii);
    phi_eta1 = testfunction_Legendre(eta1,ell_eta1);
    phi_eta2 = testfunction_Legendre(eta2,ell_eta2);
    zeta = phi_xii.*phi_eta1.*phi_eta2;
end

end

