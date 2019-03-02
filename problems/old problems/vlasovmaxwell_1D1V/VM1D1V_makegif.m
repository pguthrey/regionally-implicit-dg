function [data] = VM1D1V_makegif(fcoeffs,E1coeffs,tnow,data,nstep)
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


data.VMplotname{1} = 'Electron Density';
data.VMplotname{2} = 'Electric Field 1';

if nstep == 0
	delay = 1;
elseif tnow == data.VMFinalTime
	delay = 1;
else 
	delay = 0.1;
end

    disp(['making image , nstep is ' num2str(nstep)])

%if nstep == 1

foldername = [data.problemname '_' num2str(data.Nv1) '_M' num2str(data.M) '_r' num2str(data.r_param)];
mkdir('results/',foldername)
fileattrib(['results/' foldername],'+w')
path = ['results/' foldername '/'];
for eqn = 1
    varpath1{eqn} = ['gifs/' data.VMplotname{eqn} '_view1/' ];
    varpath2{eqn} = ['gifs/' data.VMplotname{eqn} '_topdown/'];
    varpath3{eqn} = ['gifs/' data.VMplotname{eqn} '_contour/'];
    varpath4{eqn} = ['stills/' data.VMplotname{eqn} '_view1/'];
    varpath5{eqn} = ['stills/' data.VMplotname{eqn} '_topdown/'];
    varpath6{eqn} = ['stills/' data.VMplotname{eqn} '_contour/'];

    mkdir(path,varpath1{eqn})
    mkdir(path,varpath2{eqn})
    mkdir(path,varpath3{eqn})
    mkdir(path,varpath4{eqn})
    mkdir(path,varpath5{eqn})
    mkdir(path,varpath6{eqn})

    gifname1{eqn} = [path varpath1{eqn} data.VMplotname{eqn} '.png'];
    gifname2{eqn} = [path varpath2{eqn} data.VMplotname{eqn} '.png'];
    gifname3{eqn} = [path varpath3{eqn} data.VMplotname{eqn} '.png'];
    stillname1{eqn} = [path varpath4{eqn} data.VMplotname{eqn} '_' num2str(nstep) '.png'];
    stillname2{eqn} = [path varpath5{eqn} data.VMplotname{eqn} '_' num2str(nstep) '.png'];
    stillname3{eqn} = [path varpath6{eqn} data.VMplotname{eqn} '_' num2str(nstep) '.png'];

end

for eqn = 2
    varpath1{eqn} = ['gifs/' data.VMplotname{eqn} '/' ];
    varpath4{eqn} = ['stills/' data.VMplotname{eqn} '/'];
%%{
    mkdir(path,varpath1{eqn})
   mkdir(path,varpath4{eqn})
%}
    gifname1{eqn} = [path varpath1{eqn} data.VMplotname{eqn} '.png'];
    stillname1{eqn} = [path varpath4{eqn} data.VMplotname{eqn} '_' num2str(nstep) '.png'];
end

%end


fcoeffs_1D1v = fcoeffs;

xplot = zeros(data.Nx*length(data.vmplotlocs),data.Nv1*length(data.vmplotlocs));
v1plot = zeros(data.Nx*length(data.vmplotlocs),data.Nv1*length(data.vmplotlocs));
DGplot = zeros(data.Nplotvars,data.Nx*length(data.vmplotlocs),data.Nv1*length(data.vmplotlocs));

for ix = 1:data.Nx
for iv1 = 1:data.Nv1
    xvals = data.xcenters(ix)+data.deltax*data.vmplotlocs/2;
    v1vals = data.v1centers(iv1)+data.deltav1*data.vmplotlocs/2;
    Ixii = (ix-1)*length(data.vmplotlocs);
    Ieta = (iv1-1)*length(data.vmplotlocs);
    for ixii = 1:length(data.vmplotlocs)
    for ieta = 1:length(data.vmplotlocs)
        xplot(Ixii+ixii,Ieta+ieta) = xvals(ixii);
        v1plot(Ixii+ixii,Ieta+ieta) = v1vals(ieta);
        temp = data.zeta1D1Vplot(:,:,ixii,ieta)*fcoeffs_1D1v(:,ix,iv1);
        DGplot(:,Ixii+ixii,Ieta+ieta) = temp;
    end
    end    
end
end
    
eqn = 1;

    figure(eqn)   
    set(gcf,'color','w')   
    zplot(:,:) = DGplot(eqn,:,:);
    
    clf 
    hold off
    vmax = data.boundsplot(eqn,2);
    vmin = data.boundsplot(eqn,1);
    v = linspace(vmin,vmax,100);
    contour(xplot,v1plot,zplot,v)
    axis([data.x_lb data.x_ub data.boundsv1v2v3(1,1:2)])
    colorbar        
    title([data.VMplotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    %shading interp
    xlabel(['x_1 N = ' num2str(data.Nx)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv1)],'FontSize',font)
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
    %}
    print(stillname3{eqn},'-dpng')      
    
    clf 
    hold off
    surf(xplot,v1plot,zplot)
    view(3)
    axis([data.x_lb data.x_ub data.v1_lb data.v1_ub data.vmboundsplot(eqn,:) data.vmboundsplot(eqn,:)])
    colorbar
    title([data.VMplotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nx)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv1)],'FontSize',font)
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


clear DGplot xplot
for eqn = 2
    DGsoln = E1coeffs;


    for ix = 1:data.Nx
        xvals = data.xcenters(ix)+data.deltax*data.vmplotlocs/2;
        Ixii = (ix-1)*length(data.vmplotlocs);
        for ixii = 1:length(data.vmplotlocs)
            xplot(Ixii+ixii,1) = xvals(ixii);
            temp = data.zeta1Dplot(1,:,ixii)*DGsoln(:,ix);
            DGplot(1,Ixii+ixii) = temp;
        end
    end
    

    figure(eqn)
    clf 
    set(gcf,'color','w')
    hold off
    plot(xplot,DGplot)
    axis([data.x_lb data.x_ub data.vmboundsplot(eqn,:)])
    title([data.VMplotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
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