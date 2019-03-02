function [data] = plotting_makegif_1D(DGsolution,tnow,data,nstep)
% written by Pierson Guthrey

font = 16;


if nstep == 0
	delay = 1;
elseif tnow == data.Tfinal
	delay = 1;
else 
	delay = 0.1;
end


foldername = [data.problemname '_M' num2str(data.M) '_r' num2str(data.r_param) '_' num2str(data.Nv1)];

mkdir('results/',foldername)
fileattrib(['results/' foldername],'+w')
path = ['results/' foldername '/'];
for eqn = 1:data.Nplotvars 
    %data.gifname{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved.gif'];
    %data.gifname2{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved_view2.gif'];

    varpath1{eqn} = ['gifs/' data.plotname{eqn} '_view1/' ];
    varpath4{eqn} = ['stills/' data.plotname{eqn} '_view1/'];

    mkdir(path,varpath1{eqn})
    mkdir(path,varpath4{eqn})

    gifname1{eqn} = [path varpath1{eqn} data.plotname{eqn} '.png'];
    stillname1{eqn} = [path varpath4{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];
end

%end

DGsnapshot = DGsolution;

v1plot = zeros(data.Nv1*length(data.plotlocs),1);
DGplot = zeros(data.Nplotvars,data.Nv1*length(data.plotlocs),1);

for iv1 = 1:data.Nv1
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    Ixii = (iv1-1)*length(data.plotlocs);
    for ixii = 1:length(data.plotlocs)
        v1plot(Ixii+ixii,1) = v1vals(ixii);
        temp = data.vectphiplot(:,:,ixii)*DGsnapshot(:,iv1);
        DGplot(:,Ixii+ixii) = problem_cons2plot(temp,data);
    end
end

for eqn = 1:data.Nplotvars
    figure(eqn)
    clf 
    set(gcf,'color','w')
    zplot(:,:) = DGplot(eqn,:,:);
    hold off
    plot(v1plot,zplot)
    colorbar        
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
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
       
    
end






end