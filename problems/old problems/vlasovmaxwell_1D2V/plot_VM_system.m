function [  ] = plot_VM_system(fcoeffs,E1coeffs,A1coeffs,A2coeffs,n,data)
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
for eqn = 1:4
    figure%('OuterPosition',[50 50 900 900]);
    set(gcf,'Color','w')
    clf    
    hold on
end


fnaut = @(x,v1,v2) DGeval_1D2V(fcoeffs,x,v1,v2,data);
[ f1Dcoeffs ] = project_1D2V_to_1D(fnaut,data);

E2coeffs = A1coeffs + A2coeffs;
B3coeffs = A1coeffs - A2coeffs;    
Nx = data.Nx;

for ix = 1:Nx

    xvals = data.xcenters(ix)+data.deltax*data.plotlocs/2;
    xplot = xvals;

    for ixii = 1:data.P        
        fplot(:,ixii) = data.zeta1Dplot(1,:,ixii)*f1Dcoeffs(:,ix);
        E1plot(:,ixii) = data.zeta1Dplot(1,:,ixii)*E1coeffs(:,ix,n);
        E2plot(:,ixii) = data.zeta1Dplot(1,:,ixii)*E2coeffs(:,ix,n);
        B3plot(:,ixii) = data.zeta1Dplot(1,:,ixii)*B3coeffs(:,ix,n);
    end
    
    figure(1)
    plot(xplot,fplot,'-b')
    %title(['t = ' num2str(time)],'FontSize',font)
    %axis([boundsx boundszc(eqn,:)])
    hold on
    figure(2)
    plot(xplot,E1plot,'-b')
    hold on
    figure(3)
    plot(xplot,E2plot,'-b')
    hold on
    figure(4)
    plot(xplot,B3plot,'-b')
    hold on

end

    figure(1)
    title('f')
    figure(2)
    title('E1')
    figure(3)
    title('E2')
    figure(4)
    title('B3')    

for eqn = 1:4
    figure(eqn)
    xlabel(['x   Nx = ' num2str(Nx)],'FontSize',font)
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
end
drawnow



