function [data] = initialize_mesh(data)
% Initialize mesh parameters
% written by Pierson Guthrey

data.v1_lb = data.boundsv1v2v3(1);
data.v1_ub = data.boundsv1v2v3(2);
data.v2_lb = data.boundsv1v2v3(3);
data.v2_ub = data.boundsv1v2v3(4);
data.v3_lb = data.boundsv1v2v3(5);
data.v3_ub = data.boundsv1v2v3(6);

%Derived quanities
v1endpts = linspace(data.v1_lb,data.v1_ub,data.Nv1+1);
data.v1centers = (v1endpts(1:end-1)+v1endpts(2:end))/2;
data.deltav1 = (data.v1_ub-data.v1_lb)/data.Nv1;
v2endpts = linspace(data.v2_lb,data.v2_ub,data.Nv2+1);
data.v2centers = (v2endpts(1:end-1)+v2endpts(2:end))/2;
data.deltav2 = (data.v2_ub-data.v2_lb)/data.Nv2;
v3endpts = linspace(data.v3_lb,data.v3_ub,data.Nv3+1);
data.v3centers = (v3endpts(1:end-1)+v3endpts(2:end))/2;
data.deltav3 = (data.v3_ub-data.v3_lb)/data.Nv3;
deltat1 = data.deltav1/data.Fspeedmax*data.cfl;
deltat2 = data.deltav2/data.Gspeedmax*data.cfl;
deltat3 = data.deltav3/data.Hspeedmax*data.cfl;
deltat = min([deltat1 deltat2 deltat3]);
data.nuv1 = deltat/data.deltav1;
data.nuv2 = deltat/data.deltav2;
data.nuv3 = deltat/data.deltav3;
end    












       
        
