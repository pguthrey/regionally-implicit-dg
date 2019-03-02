function [data] = plotting_makegif_2Dto1D(DGsolution,tnow,data,nstep)
%  Plots the DG solution at the current time
% written by Pierson Guthrey

font = 16;


if nstep == 0
	delay = 1;
elseif tnow == data.Tfinal
	delay = 1;
else 
	delay = 0.1;
end

  %  disp(['making image , nstep is ' num2str(nstep)])


if nstep == 0
   for eqn = 1:data.Nplotvars 
         %fig = figure('OuterPosition',[50 50 900 900]);
        fig = figure;
        data.fignum{eqn} = fig.Number;
   end
end

foldername = [data.problemname '_' num2str(data.Nv1) '_M' num2str(data.M) '_r' num2str(data.r_param)];

mkdir('results/',foldername)
fileattrib(['results/' foldername],'+w')
path = ['results/' foldername '/'];
for eqn = 1:data.Nplotvars 
    varpath1{eqn} = ['gifs_1D/' data.plotname{eqn} '/'];
    varpath2{eqn} = ['stills_1D/' data.plotname{eqn} '/'];
    mkdir(path,varpath1{eqn})
    mkdir(path,varpath2{eqn})
    data.gifname{eqn} = [path varpath1{eqn} data.plotname{eqn} '.png'];
    data.stillname{eqn} = [path varpath2{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];
    %fig = figure('OuterPosition',[50 50 900 900]);
    set(gcf,'Color','w')
end

%end

ivhalf = ceil(data.Nv2/2);

DGsnapshot(:,:) = DGsolution(:,:,ivhalf);
    
v1plot = zeros(data.Nv1*length(data.plotlocs),1);
DGplot = zeros(data.Nplotvars,data.Nv1*length(data.plotlocs));

for iv1 = 1:data.Nv1
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    Ixii = (iv1-1)*length(data.plotlocs);
    for ixii = 1:length(data.plotlocs)
        v1plot(Ixii+ixii,1) = v1vals(ixii);
        temp = data.vectphi_2Dto1D(:,:,ixii)*DGsnapshot(:,iv1);        
        DGplot(:,Ixii+ixii) = problem_cons2plot(temp,data);
    end    
end

for eqn = 1:data.Nplotvars
    figure(data.Nplotvars+eqn)
    clf 
    set(gcf,'color','w')
    zplot = DGplot(eqn,:);
    plot(v1plot,zplot)
    axis([data.boundsv1v2v3(1,1:2) data.boundsplot(eqn,:)])
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    set(gca,'FontSize',font)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if nstep == 0
        imwrite(imind,cm,data.gifname{eqn},'gif','Loopcount',inf,'DelayTime',delay);
    else%if ~mod(nstep,frameskip)
        imwrite(imind,cm,data.gifname{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    print(data.stillname{eqn},'-dpng')
           
end
