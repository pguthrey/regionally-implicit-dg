function [ DGcoeffs_new ] = project_1D2V_to_1D2V(fxv1v2,data)

%h1 = waitbar(0,'Advecting 1D2V to 1D2V','OuterPosition', [400 100 300 75]);
%{
data.locs = data.locs;
P = data.P;
Nx = data.Nx;
Nv1 = data.Nv1;
Nv2 = data.Nv2;
locs = data.locs;
xcenters = data.xcenters;
v1centers = data.v1centers;
v2centers = data.v2centers;
deltax = data.deltax;
deltav1 = data.deltav1;
deltav2 = data.deltav2;
theta1D2V = data.theta1D2V;
zeta1D2V = data.zeta1D2V;
%}

%%{
normphi = 8;
DGcoeffs_new = zeros(data.theta1D2V,data.Nx,data.Nv1,data.Nv2);
for ix = 1:data.Nx
for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    temp = zeros(data.theta1D2V,1);
    for xquad = 1:data.P
    for v1quad = 1:data.P
    for v2quad = 1:data.P
        xii = data.locs(xquad);
        eta1 = data.locs(v1quad);
        eta2 = data.locs(v2quad);
        x = data.xcenters(ix)+xii*data.deltax/2;
        v1 = data.v1centers(iv1)+eta1*data.deltav1/2;
        v2 = data.v2centers(iv2)+eta2*data.deltav2/2;
        ftilde = fxv1v2(x,v1,v2);       
        zetatrans = data.zeta1D2V(1,:,xquad,v1quad,v2quad)';
        wgt = data.wgts3D(xquad,v1quad,v2quad);
        temp = temp + zetatrans*ftilde*wgt/normphi;
    end
    end
    end
    %keyboard
    DGcoeffs_new(:,ix,iv1,iv2) = temp;
end
end
%waitbar(ix/Nx,h1)
end
%}

%{
M = 400
[xgrid,v1grid] = meshgrid(linspace(-6+10*eps,6,M));
for ix = 1:M
for iv = 1:M    
    x = xgrid(ix,iv);
    v1 = v1grid(ix,iv);
    DG(ix,iv) = fxv1v2(x,v1,0); 
end
end
close all
surf(DG)
shading interp
%}
%A(:,:) = DGcoeffs_new(1,:,:,51); close all ; surf(A) ; shading interp

%{
DGcoeffs_new = zeros(data.theta1D2V,data.Nx,data.Nv1,data.Nv2);
for ix = 1:data.Nx
for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
        x = data.xcenters(ix);
        v1 = data.v1centers(iv1);
        ftilde = fxv1v2(x,v1,0);
        DGcoeffs_new(:,ix,iv1,iv2) = ftilde;
end
end
%waitbar(ix/Nx,h1)
end

%}

