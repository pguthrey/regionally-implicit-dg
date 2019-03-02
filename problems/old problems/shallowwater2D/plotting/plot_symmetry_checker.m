function [filename1] = plot_symmetry_checker(DGsnapshot,time,timename,path,filename,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

font = 16;

%load('data/data_parameters')
%load('data/data_plotters')
%load('data/data_quadrature')

fig1 = figure('OuterPosition',[0 0 100+800*data.Neqns 900]);
data.fignum = fig1.Number;
figure(data.fignum)
set(gcf,'Color','w')
clf    
hold on
colorbar

filename1 =[path filename '_T' timename '_symmetry_check.gif'];

r = zeros(data.Neqns,1);
for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    v2vals = data.v2centers(iv2)+data.deltav2*data.plotlocs/2;
    for ixii = 1:data.P
    for ieta = 1:data.P
        v1plot = v1vals(ixii);
        v2plot = v2vals(ieta);
        temp = data.vectphiplot(:,:,ixii,ieta)*DGsnapshot(:,iv1,iv2);
        r = sqrt(v1plot^2+v2plot^2);
        for eqn = 1:data.Neqns
            subplot(1,data.Neqns,eqn)
            plot(r,temp(eqn),'.');
        end            
    end
    end
    view(2)
    xmax = max([abs(boundsv1v2(1)) abs(boundsv1v2(2))]);
    ymax = max([abs(boundsv1v2(3)) abs(boundsv1v2(4))]);
    rmax = sqrt(xmax^2+ymax^2);
    %axis([0 rmax boundsq(eqn,1:2)])
    hold on
    drawnow
end
end

subplot(1,data.Neqns,1)
title(['t = ' num2str(time)],'FontSize',font)

for eqn = 1:data.Neqns
    subplot(1,data.Neqns,eqn)
    shading interp
    xlabel(['Distance of sample point ||x|| '],'FontSize',font)
    ylabel(['Solution values'],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
end

drawnow
frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

imwrite(imind,cm,filename1,'gif');


end

