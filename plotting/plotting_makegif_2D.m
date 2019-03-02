function [data] = plotting_makegif_2D(DGsolution,tnow,data,nstep)
% written by Pierson Guthrey

font = 16;
%{
deltatvals = diff(data.Time);
if all(deltatvals > .06)
    delay = [2 deltatvals 2 2];
    frameskip=1;
else 
    delay = [2 deltatvals*0.06/min(deltatvals)*3 2 2];    
    frameskip=1;
end
%}

if nstep == 0
	delay = 1;
elseif tnow == data.Tfinal
	delay = 1;
else 
	delay = 0.1;
end

   % disp(['making image , nstep is ' num2str(nstep)])

%if nstep == 1

foldername = [data.problemname '_M' num2str(data.M) '_r' num2str(data.r_param) '_' num2str(data.Nv1)];

mkdir('results/',foldername)
fileattrib(['results/' foldername],'+w')
path = ['results/' foldername '/'];
for eqn = 1:data.Nplotvars 
    %data.gifname{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved.gif'];
    %data.gifname2{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved_view2.gif'];

    varpath1{eqn} = ['gifs/' data.plotname{eqn} '_view1/' ];
    varpath2{eqn} = ['gifs/' data.plotname{eqn} '_topdown/'];
    varpath3{eqn} = ['gifs/' data.plotname{eqn} '_contour/'];
    varpath4{eqn} = ['stills/' data.plotname{eqn} '_view1/'];
    varpath5{eqn} = ['stills/' data.plotname{eqn} '_topdown/'];
    varpath6{eqn} = ['stills/' data.plotname{eqn} '_contour/'];

    mkdir(path,varpath1{eqn})
    mkdir(path,varpath2{eqn})
   mkdir(path,varpath3{eqn})
   mkdir(path,varpath4{eqn})
   mkdir(path,varpath5{eqn})
   mkdir(path,varpath6{eqn})

    gifname1{eqn} = [path varpath1{eqn} data.plotname{eqn} '.png'];
    gifname2{eqn} = [path varpath2{eqn} data.plotname{eqn} '.png'];
    gifname3{eqn} = [path varpath3{eqn} data.plotname{eqn} '.png'];
    stillname1{eqn} = [path varpath4{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];
    stillname2{eqn} = [path varpath5{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];
    stillname3{eqn} = [path varpath6{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];

end

%end

DGsnapshot(:,:,:) = DGsolution(:,:,:);

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
    view(3)
    %axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
    colorbar        
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %%{
    if nstep == 0
        imwrite(imind,cm,gifname1{eqn},'gif','Loopcount',inf,'DelayTime',delay);
    else%if ~mod(nstep,frameskip)
        imwrite(imind,cm,gifname1{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    %}
    print(stillname1{eqn},'-dpng')
       
    view(2)
    %axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if nstep == 0
        imwrite(imind,cm,gifname2{eqn},'gif','Loopcount',inf,'DelayTime',1);
    else%if ~mod(n,frameskip)
        imwrite(imind,cm,gifname2{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    print(stillname2{eqn},'-dpng')

    %{
    clf 
    hold off
    vmax = data.boundsplot(eqn,2);
    vmin = data.boundsplot(eqn,1);
    v = linspace(vmin,vmax,100);
    contour(v1plot,v2plot,zplot,v)
    %view(2)
    %axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
    colorbar        
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    %shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);    
    
    %%{
    if nstep == 0
        imwrite(imind,cm,gifname3{eqn},'gif','Loopcount',inf,'DelayTime',1);
    else%if ~mod(n,frameskip)
        imwrite(imind,cm,gifname3{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    print(stillname3{eqn},'-dpng')
    %}
    
end

end
