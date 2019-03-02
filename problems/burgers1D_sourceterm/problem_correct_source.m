function [correct_cell,data] = problem_correct_source(correct_cell,qstar_1D,auxiliary,data,cellcenter)

vectpsi = data.vectpsi;

vectphi_Trans = data.vectphi_Trans;

Pdp1 = data.Pdp1;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
beta = data.beta;
nuv1 = data.nuv1;

for k = 1:Pdp1
    tquad = Dp1list(k,1);
    v1quad = Dp1list(k,2);
    weight = Dp1quadwgts(k);
    
    wstar = vectpsi(:,:,tquad,v1quad)*qstar_1D;
        
    correct_cell = correct_cell - (nuv1*weight)*vectphi_Trans(:,:,v1quad)*beta*wstar;
end


end

