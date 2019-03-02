function [  ] = plot_2D_system(DGsnapshot,time,data)
% data.Plots the DG solution at the current time
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    x,y        : cell centers (these should be %loaded)
%           DGpast       : coefficients from the previous timestep
% OUTdata.PUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------
font = 16;

for eqn = 1:data.Neqns
    figure(1)
    clf
end

v1plot = zeros(data.Pplot,data.Pplot);
v2plot = zeros(data.Pplot,data.Pplot);
DGplot = zeros(data.Neqns,data.Pplot,data.Pplot);
xii_full = zeros(data.Pplot*data.Nv1,data.Pplot*data.Nv2);
eta_full = zeros(data.Pplot*data.Nv1,data.Pplot*data.Nv2);
DGplot_full = zeros(data.Neqns,data.Pplot*data.Nv1,data.Pplot*data.Nv2);

for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2        
    v1vals = data.v1centers(iv1)+data.deltav1*data.Dplotlocs(:,1)/2;
    v2vals = data.v2centers(iv2)+data.deltav2*data.Dplotlocs(:,2)/2;
        
    for vplot = 1:data.PDplot
        v1celli = data.plotlist(vplot,1);
        v2celli = data.plotlist(vplot,2);
        v1plot(v1celli,v2celli) = v1vals(vplot);
        v2plot(v1celli,v2celli) = v2vals(vplot);
        temp = data.vectphiplot(:,:,v1celli,v2celli)*DGsnapshot(:,iv1,iv2);
        DGplot(:,v1celli,v2celli) = temp;
    end    
    %keyboard

    thing1= (1:data.Pplot)+(iv1-1)*data.Pplot;
    thing2= (1:data.Pplot)+(iv2-1)*data.Pplot;
    
    xii_full(thing1,thing2) = v1plot;
    eta_full(thing1,thing2) = v2plot;
    DGplot_full(:,thing1,thing2) = DGplot;
end
end

for eqn = 1:data.Neqns
    figure(1)
    qplot(:,:) = DGplot_full(eqn,:,:);
    hold on
    surf(xii_full,eta_full,qplot)		 
    shading interp
end
    drawnow
    %axis([data.boundsv1v2v3 data.boundsq(eqn,:)])

for eqn = 1:data.Neqns
    %figure(data.fignum{eqn})
    figure(1)
    title([data.varname{eqn} ' , t = ' num2str(time)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    colorbar
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    drawnow
    view(2)
end
pause(2)
end