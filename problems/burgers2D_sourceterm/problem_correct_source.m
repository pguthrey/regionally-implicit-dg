function [correct_cell,auxiliary,data] = problem_correct_source(correct_cell,qstar,auxiliary,data,cellcenter,cellindex)

vectphi_Trans = data.vectphi_Trans;
vectpsi = data.vectpsi;

%%{
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
%}

vectphi = data.vectphi;
Pd = data.Pd;
space_dims = data.space_dims;
Dlist = data.Dlist;
Dquadwgts = data.Dquadwgts;


beta = data.beta;

%qstar_phi = 2*qstar(data.BO(:,1)==1,1);

v2quad = 1;
v3quad = 1;
%{
for k = 1:Pd
        v1quad = Dlist(k,1);
        if space_dims >= 2
            v2quad = Dlist(k,2);
            if space_dims >= 3
                v3quad = Dlist(k,3);
            end
        end
        weight = Dquadwgts(k);
        wstar = vectpsi(:,:,v1quad,v2quad,v3quad)*qstar_phi;
        source = beta*wstar;
        correct_cell = correct_cell - vectphi_Trans(:,:,v1quad,v2quad,v3quad)*source*weight;
end
%}

%%{
for k = 1:Pdp1
        tquad = Dp1list(k,1);
        v1quad = Dp1list(k,2);
        v2quad = Dp1list(k,3);
        if space_dims >= 3
            v3quad = Dp1list(k,4);
        else
            v3quad = 1;
        end
        weight = Dp1quadwgts(k);
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        source = beta*wstar;
        correct_cell = correct_cell - vectphi_Trans(:,:,v1quad,v2quad,v3quad)*source*weight;
end
%}
end

