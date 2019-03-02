function [ DGcoeffs_new ] = project_1D2V_to_1D(fxv1v2,data)
%load('vlasovmaxwell/data.appdata')
%h1 = waitbar(0,'Advecting 1D2V to 1D','OuterPosition', [400 100 300 75]);

data.locs = data.locs;
P = data.P;
Nx = data.Nx;
locs = data.locs;
xcenters = data.xcenters;
v1centers = data.v1centers;
v2centers = data.v2centers;
deltax = data.deltax;
deltav1 = data.deltav1;
deltav2 = data.deltav2;
zeta1D = data.zeta1D;
wgts3D = data.wgts3D;
DGcoeffs_new = zeros(data.theta1D,Nx);
normphi = 2;
for ix = 1:Nx
    temp = zeros(data.theta1D,1);
    for iv1 = 1:data.Nv1
    for iv2 = 1:data.Nv2
        for v1quad = 1:P
        for v2quad = 1:P
        for xquad = 1:P
            xii = locs(xquad);
            eta1 = locs(v1quad);
            eta2 = locs(v2quad);
            x = xcenters(ix)+xii*deltax/2;
            v1 = v1centers(iv1)+eta1*deltav1/2;
            v2 = v2centers(iv2)+eta2*deltav2/2;
            zetatrans = zeta1D(1,:,xquad)';
            ftilde = fxv1v2(x,v1,v2);
            wgt = wgts3D(xquad,v1quad,v2quad);
            temp = temp + zetatrans*ftilde*wgt/normphi;
        end
        end
        end
    end
    end
    DGcoeffs_new(:,ix) = temp;
    %waitbar(ix/Nx,h1)
end
%delete(h1)

end
