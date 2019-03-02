function [bigP] = compute_viscocity(DGsolution_old,data)


n_averages(:,:,:)  = DGsolution_old(0*data.Ls+1,:,:,:);
m1_averages(:,:,:) = DGsolution_old(1*data.Ls+1,:,:,:);
m2_averages(:,:,:) = DGsolution_old(2*data.Ls+1,:,:,:);
m3_averages(:,:,:) = DGsolution_old(3*data.Ls+1,:,:,:);

u1_averages = m1_averages./n_averages;
u2_averages = m2_averages./n_averages;
u3_averages = m3_averages./n_averages;

dx = data.deltav1;
dy = data.deltav2;
dz = data.deltav3;

u1x = (circshift(u1_averages,-1,1) - circshift(u1_averages,1,1))/(2*dx); 
u1y = (circshift(u1_averages,-1,2) - circshift(u1_averages,1,2))/(2*dy); 
u1z = (circshift(u1_averages,-1,3) - circshift(u1_averages,1,3))/(2*dz);

u2x = (circshift(u2_averages,-1,1) - circshift(u2_averages,1,1))/(2*dx); 
u2y = (circshift(u2_averages,-1,2) - circshift(u2_averages,1,2))/(2*dy); 
u2z = (circshift(u2_averages,-1,3) - circshift(u2_averages,1,3))/(2*dz);

u3x = (circshift(u3_averages,-1,1) - circshift(u3_averages,1,1))/(2*dx); 
u3y = (circshift(u3_averages,-1,2) - circshift(u3_averages,1,2))/(2*dy); 
u3z = (circshift(u3_averages,-1,3) - circshift(u3_averages,1,3))/(2*dz);

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

for iv1 = 1:Nv1 
for iv2 = 1:Nv2 
for iv3 = 1:Nv3 
    
    n0 = n_averages(iv1,iv2,iv3);
    [params] = cell_average_to_params(n0,data);
    eta = params.eta;
    xi = params.xi;
    
    %first row
    ugrads(1,1,:,:,:) = eta*(u1x + u1x - 2/3*u1x ) + xi*u1x ;
    ugrads(1,2,:,:,:) = eta*(u1y + u2x) ;
    ugrads(1,3,:,:,:) = eta*(u1z + u3x) ;
    %second row
    ugrads(2,1,:,:,:) = eta*(u2x + u1y) ;
    ugrads(2,2,:,:,:) = eta*(u2y + u2y - 2/3*u2y ) + xi*u2y ;
    ugrads(2,3,:,:,:) = eta*(u2z + u3y) ;
    %third row
    ugrads(3,1,:,:,:) = eta*(u3x + u1z) ;
    ugrads(3,2,:,:,:) = eta*(u3y + u2z) ;
    ugrads(3,3,:,:,:) = eta*(u3z + u3z - 2/3*u3z ) + xi*u3z ;
            
end
end
end

ugrads1(:,:,:,:) = ugrads(:,1,:,:,:);
ugrads1_dx = (circshift(ugrads1,-1,2) - circshift(ugrads1,1,2))/(2*dx); 
ugrads2(:,:,:,:) = ugrads(:,2,:,:,:);
ugrads2_dy = (circshift(ugrads2,-1,3) - circshift(ugrads2,1,3))/(2*dy); 
ugrads3(:,:,:,:) = ugrads(:,3,:,:,:);
ugrads3_dz = (circshift(ugrads3,-1,4) - circshift(ugrads3,1,4))/(2*dz); 

bigP = ugrads1_dx + ugrads2_dy + ugrads3_dz;

end

