function [data] = problem_makegif(DGsolution,tnow,data,nstep)
% data.Plots the DG solution at the current time
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y        : cell centers (these should be %loaded)
%           DGpast       : coefficients from the previous timestep
% OUTdata.PUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
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
delay = .1;

    disp(['making image , nstep is ' num2str(nstep)])

%if nstep == 1

foldername = [data.problemname '_' num2str(data.Nv1) '_M' num2str(data.M) '_r' num2str(data.r_param)];

mkdir('results/',foldername)
fileattrib(['results/' foldername],'+w')
path = ['results/' foldername '/'];
for eqn = 1:data.Nplotvars 
    %data.gifname{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved.gif'];
    %data.gifname2{end+1} =[data.problemname '_' data.plotname{eqn} '_conserved_view2.gif'];

    varpath1{eqn} = ['gifs/' data.plotname{eqn} '_view1/' ];
    varpath2{eqn} = ['gifs/' data.plotname{eqn} '_topdown/'];
    varpath3{eqn} = ['stills/' data.plotname{eqn} '_view1/'];
    varpath4{eqn} = ['stills/' data.plotname{eqn} '_topdown/'];

    mkdir(path,varpath1{eqn})
    mkdir(path,varpath2{eqn})
   mkdir(path,varpath3{eqn})
   mkdir(path,varpath4{eqn})

    gifname1{eqn} = [path varpath1{eqn} data.plotname{eqn} '.png'];
    gifname2{eqn} = [path varpath2{eqn} data.plotname{eqn} '.png'];
    stillname1{eqn} = [path varpath3{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];
    stillname2{eqn} = [path varpath4{eqn} data.plotname{eqn} '_' num2str(nstep) '.png'];

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
    zplot(:,:) = DGplot(eqn,:,:);
    hold on
    surf(v1plot,v2plot,zplot)
    view(3)
    axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
    colorbar        
end
    

for eqn = 1:data.Nplotvars
    figure(eqn)
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    colorbar
    axis([data.boundsv1v2v3(1,1:4) data.boundsplot(eqn,:) data.boundsplot(eqn,:)])
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
    axis([data.boundsv1v2v3(1,1:4) data.boundsq(eqn,:) data.boundsq(eqn,:)])
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %%{
    if nstep == 0
        imwrite(imind,cm,gifname2{eqn},'gif','Loopcount',inf,'DelayTime',delay);
    else%if ~mod(n,frameskip)
        imwrite(imind,cm,gifname2{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    %}
    print(stillname2{eqn},'-dpng')
    
end
end