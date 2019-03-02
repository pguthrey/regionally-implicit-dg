function [gifname,gifname2] = plot_2D_system_makegif_conserved(DGsolution,path,filename,data)
% data.Plots the DG solution at the current time
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y        : cell centers (these should be %loaded)
%           DGpast       : coefficients from the previous timestep
% OUTdata.PUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
font = 16;
close all

deltatvals = diff(data.Time);
if all(deltatvals > .06)
    delay = [2 deltatvals 2 2];
    frameskip=1;
else 
    delay = [2 deltatvals*0.06/min(deltatvals)*3 2 2];    
    frameskip=1;
end

gifname = {};
gifname2 = {};

for eqn = 1:data.Neqns 
    gifname{end+1} =[path filename '_' data.varname{eqn} '_conserved.gif'];
    gifname2{end+1} =[path filename '_' data.varname{eqn} '_conserved_view2.gif'];
    fig = figure('OuterPosition',[50 50 900 900]);
    set(gcf,'Color','w')
    data.fignum{eqn} = fig.Number;
end

k = 0;
for n = [1 1:data.Nt data.Nt]
    k = k+1;
    DGsnapshot(:,:,:) = DGsolution(:,:,:,n);
    clf 
    figure(data.fignum{eqn})
    hold on
for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    v2vals = data.v2centers(iv2)+data.deltav2*data.plotlocs/2;
    for ixii = 1:length(data.plotlocs)
    for ieta = 1:length(data.plotlocs)
        v1plot(ixii,ieta) = v1vals(ixii);
        v2plot(ixii,ieta) = v2vals(ieta);
        temp = data.vectphiplot(:,:,ixii,ieta)*DGsnapshot(:,iv1,iv2);
        DGplot(:,ixii,ieta) = temp;
    end
    end    
    
    for eqn = 1:data.Neqns
        figure(data.fignum{eqn})
        zplot(:,:) = DGplot(eqn,:,:);
        hold on
        surf(v1plot,v2plot,zplot)
        view(3)
		%axis([boundsv1v2 boundsq(eqn,:)])
        colorbar        
    end
    
end
end

for eqn = 1:data.Neqns
    figure(data.fignum{eqn})
    title([data.varname{eqn} ' , t = ' num2str(data.Time(n))],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    colorbar
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
%{
    if n == 1
        imwrite(imind,cm,gifname{eqn},'gif','Loopcount',inf,'DelayTime',delay(1));
    elseif ~mod(n,frameskip)
        imwrite(imind,cm,gifname{eqn},'gif','WriteMode','append','DelayTime',delay(k));
    end
    %}
    imwrite(imind,cm,[gifname{eqn} num2str(step)],'gif','Loopcount',inf,'DelayTime',delay(1));

    view(2)
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if n == 1
        imwrite(imind,cm,gifname2{eqn},'gif','Loopcount',inf,'DelayTime',delay(1));
    elseif ~mod(n,frameskip)
        imwrite(imind,cm,gifname2{eqn},'gif','WriteMode','append','DelayTime',delay(k));
    end
end
end