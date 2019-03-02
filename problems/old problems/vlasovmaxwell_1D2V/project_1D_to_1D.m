function [ DGcoeffs_new ] = project_1D_to_1D(fx,data)
%load('vlasovmaxwell/data.appdata')

%h1 = waitbar(0,'Advecting 1D to 1D','OuterPosition', [600 100 300 75]);

data.locs = data.locs;
P = data.P;
Nx = data.Nx;
locs = data.locs;
xcenters = data.xcenters;
deltax = data.deltax;
wgts1D = data.wgts1D;
zeta1D = data.zeta1D;
theta1D = data.theta1D;

DGcoeffs_new = zeros(theta1D,Nx);
normphi = 2;

for ix = 1:Nx
    temp = zeros(theta1D,1);
    for xquad = 1:P    
        xii = locs(xquad);
        x = xcenters(ix)+xii*deltax/2;
        zetatrans = zeta1D(1,:,xquad)';
        temp = temp + zetatrans*fx(x)*wgts1D(xquad)./normphi;
    end
    DGcoeffs_new(:,ix) = temp;
    %waitbar(ix/Nx,h1)
end


