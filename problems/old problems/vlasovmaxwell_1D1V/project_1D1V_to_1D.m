function [ DGcoeffs_new ] = project_1D1V_to_1D(fxv1v2,data)
%load('vlasovmaxwell/data.appdata')
%h1 = waitbar(0,'Advecting 1D2V to 1D','OuterPosition', [400 100 300 75]);


DGcoeffs_new = zeros(data.theta1D,data.Nx);
normphi = 2;
for ix = 1:data.Nx
    temp = zeros(data.theta1D,1);
    for iv1 = 1:data.Nv1
        for v1quad = 1:data.P
        for xquad = 1:data.P
            xii = data.locs(xquad);
            eta1 = data.locs(v1quad);
            x = data.xcenters(ix)+xii*data.deltax/2;
            v1 = data.v1centers(iv1)+eta1*data.deltav1/2;
            zetatrans = data.zeta1D(1,:,xquad)';
            ftilde = fxv1v2(x,v1);
            wgt = data.wgts2D(xquad,v1quad);
            temp = temp + zetatrans*ftilde*wgt/normphi;
        end
        end
    end
    DGcoeffs_new(:,ix) = temp;
    %waitbar(ix/Nx,h1)
end
%delete(h1)

end
