function [data] = plotting_makegif_1D(DGsolution,... %
                                        tnow,... %
                                        data,... %
                                        nstep) %
% written by Pierson Guthrey

font = 16;


if nstep == 0
	delay = 1;
elseif tnow == data.Tfinal
	delay = 1;
else 
	delay = 0.1;
end

%    disp(['making image , nstep is ' num2str(nstep)])


if nstep == 0
   for eqn = 1:data.Nplotvars 
         %fig = figure('OuterPosition',[50 50 900 900]);
        fig = figure;
        data.fignum{eqn} = fig.Number;
   end
end

for eqn = 1:data.Nplotvars 
    varpath = [data.plotname{eqn} '/'];
    mkdir(path,varpath)
    data.gifname{eqn} = [path varpath data.plotname{eqn} '.png'];
    data.stillname{eqn} = [path varpath data.plotname{eqn} '_' num2str(nstep) '.png'];
    %fig = figure('OuterPosition',[50 50 900 900]);
    set(gcf,'Color','w')
end

%end

DGsnapshot(:,:,:) = DGsolution(:,:,:);
    
v1plot = zeros(data.Nv1*length(data.plotlocs),1);
DGplot = zeros(data.Nplotvars,data.Nv1*length(data.plotlocs));

for iv1 = 1:data.Nv1
    v1vals = data.v1centers(iv1)+data.deltav1*data.plotlocs/2;
    Ixii = (iv1-1)*length(data.plotlocs);
    for ixii = 1:length(data.plotlocs)
        v1plot(Ixii+ixii,1) = v1vals(ixii);    
        if data.flags(iv1)            
            whichaverage = min(max(floor(data.limiter_subcells*(data.plotlocs(ixii)+1)/2)+1,1),data.limiter_subcells);
            qval(:,1) = data.DG_averages(1,:,whichaverage,iv1);
            DGplot(:,Ixii+ixii) = problem_cons2plot(qval,data);
        else
            qval = data.vectphiplot(:,:,ixii)*DGsnapshot(:,iv1);
            DGplot(:,Ixii+ixii) = problem_cons2plot(qval,data);
        end    
    end
end

for eqn = 1:data.Nplotvars
    figure(eqn)
    hold off
    clf 
    set(gcf,'color','w')

    %for i = 1:length(v1plot)
    %    [ exact(i) ] = problem_solution(tnow,v1plot(i),data);    
    %end
    %plot(v1plot,exact,'k-');
    %hold on

    zplot = DGplot(eqn,:);
    plot(v1plot,zplot)
    hold on
    %{
    for k = 1:length(v1plot)
        exactsoln(k) = problem_solution(tnow,v1plot(k),data);
    end
    plot(v1plot,exactsoln,'--')
    hold off
%}
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

%{
keyboard
hold on 
plot(data.v1centers,data.flags,'o')
%}
end