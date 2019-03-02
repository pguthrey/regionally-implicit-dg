function [filename1] = plot_1D_system_makegif_conserved(solution,data)
% Plots the DG solution at the current time
% written by Pierson Guthrey
% -------------------------------------------------
% INPUTS    x,y        : cell centers (these should be %loaded)
%           DGpast       : coefficients from the previous timestep
% OUTPUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------
font = 16;

deltatvals = diff(Time);
if all(deltatvals > .06)
    delay = [2 deltatvals 2 2];
    frameskip=1;
else 
    delay = [2 deltatvals*0.06/min(deltatvals) 2 2];    
    frameskip=1;
end

%load('data/data_parameters')
%load('data/data_data.basis')
%load('data/data_mesh')
%load('data/data_plotEvals','data.vectphiplot','data.plotlocs')
%load('data/data_quadrature','data.locs','P')


fig1 = figure('OuterPosition',[50 50 100+800*data.Neqns 900]);
data.fignum = fig1.Number;

filename1 =[path filename '_conserved.gif'];

for n = [1:Nt Nt]
    DGsnapshot(:,:) = DGsolution(:,:,n);
    clf 
    figure(data.fignum)
    set(gcf,'Color','w')
    hold on
    
for iv1 = 1:data.Nv1

    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;

    for ixii = 1:data.P

        v1plot(ixii,1) = v1vals(ixii);

        temp = data.vectphiplot(:,:,ixii)*DGsnapshot(:,iv1);
        DGplot(:,ixii) = temp;
    end

    
    for eqn = 1:data.Neqns
        subplot(1,data.Neqns,eqn)
        zplot(:,1) = DGplot(eqn,:,:);
        plot(v1plot,zplot,'-b')
		axis([data.boundsv1v2 data.boundsq(eqn,:)])
        hold on
    end
end

subplot(1,data.Neqns,1)
title(['t = ' num2str(Time(n))],'FontSize',font)

for eqn = 1:data.Neqns
    subplot(1,data.Neqns,eqn)
    shading interp
    xlabel(['v1   data.Nv1 = ' num2str(data.Nv1)],'FontSize',font)
    
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
end
drawnow
frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

if n == 1
    imwrite(imind,cm,filename1,'gif','Loopcount',inf,'DelayTime',2);
elseif ~mod(n,frameskip)
    imwrite(imind,cm,filename1,'gif','WriteMode','append','DelayTime',.1);
end
end