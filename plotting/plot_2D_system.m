function [  ] = plot_2D_system(DGsnapshot,tnow,data,nstep)
% written by Pierson Guthrey

font = 16;
for eqn = 1:data.Nplotvars 
 

end

v1plot = zeros(data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));
v2plot = zeros(data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));
DGplot = zeros(data.Nplotvars,data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));

for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    v2vals = data.v2centers(iv2)+data.deltav2*data.plotlocs/2;
    Ixii = (iv1-1)*length(data.plotlocs);
    Ieta = (iv2-1)*length(data.plotlocs);
    for ixii = 1:length(data.plotlocs)
    for ieta = 1:length(data.plotlocs)
        v1plot(Ixii+ixii,Ieta+ieta) = v1vals(ixii);
        v2plot(Ixii+ixii,Ieta+ieta) = v2vals(ieta);
        temp = data.vectphiplot(:,:,ixii,ieta)*DGsnapshot(:,iv1,iv2);
        DGplot(:,Ixii+ixii,Ieta+ieta) = problem_cons2plot(temp,data);
    end
    end    
end
end

for eqn = 1:data.Nplotvars
    figure(eqn)
    clf 
    set(gcf,'color','w')
    zplot(:,:) = DGplot(eqn,:,:);
    hold off
    surf(v1plot,v2plot,zplot)
    view(2)
    axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
    colorbar        
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
   
end
end