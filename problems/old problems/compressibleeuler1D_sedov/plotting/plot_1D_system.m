function [  ] = plot_1D_system(DGsolution,tnow,data)
% Plots the DG solution at the current time
% written by Pierson Guthrey
% -------------------------------------------------
% INPUTS    x,y        : cell centers (these should be %loaded)
%           DGpast       : coefficients from the previous timestep
% OUTPUTS   DGcorrection : corrected coefficients         
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------
font = 16;

for eqn = 1:data.Nplotvars
    figure(eqn)
    clf
end

v1plot = zeros(data.Pplot,1);
DGplot = zeros(data.Nplotvars,data.Pplot);
xii_full = zeros(data.Pplot*data.Nv1,1);
DGplot_full = zeros(data.Nplotvars,data.Pplot*data.Nv1);

for iv1 = 1:data.Nv1
    v1vals = data.v1centers(iv1)+data.deltav1*data.Dplotlocs(:,1)/2;        
    for vplot = 1:length(data.plotlocs)
        v1plot(vplot,1) = v1vals(vplot);
        temp = data.vectphiplot(:,:,vplot)*DGsolution(:,iv1);
        DGplot(:,vplot) = problem_cons2plot(temp,data);
    end    
    %keyboard

    thing1= (1:data.Pplot)+(iv1-1)*data.Pplot;
    
    xii_full(thing1,1) = v1plot;
    DGplot_full(:,thing1) = DGplot;
end

for eqn = 1:data.Nplotvars
    figure(eqn)
    qplot = DGplot_full(eqn,:);
    plot(xii_full,qplot); 
    title([data.plotname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['Solution'],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    axis([data.boundsv1v2v3(1:2) data.boundsplot(eqn,:)])
end
    drawnow