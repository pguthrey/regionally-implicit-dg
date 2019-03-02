function [data] = plot_2D_makegif_conserved_steps(DGsolution,tnow,data,nstep)
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
    foldername = [data.problemname '_' num2str(data.Nv1)];

if nstep == 1
   for eqn = 1:data.Neqns 
         fig = figure('OuterPosition',[50 50 900 900]);
         data.fignum{eqn} = fig.Number;
         
        %[~,~,~] = rmdir(['results/' data.problemname]);
        mkdir('results/',foldername)
        %fileattrib(foldername,'+w')
         
   end
end

%if nstep == 1

path = ['results/' foldername '/'];
data.gifname = {};
data.gifname2 = {};
    for eqn = 1:data.Neqns 
        %data.gifname{end+1} =[data.problemname '_' data.varname{eqn} '_conserved.gif'];
        %data.gifname2{end+1} =[data.problemname '_' data.varname{eqn} '_conserved_view2.gif'];
    
        varpath = [data.varname{eqn} '_view1/' ];
        varpath2 = [data.varname{eqn} '_topdown/' ];
        
        [~,~,~] = rmdir([path,varpath]);
        mkdir(path,varpath)
        fileattrib([path,varpath],'+w')
        [~,~,~] = rmdir([path,varpath2]);
        mkdir(path,varpath2)
        fileattrib([path,varpath2],'+w')
            
        data.gifname{end+1} = [path varpath data.varname{eqn} '_' num2str(nstep) '.png'];
        data.gifname2{end+1} = [path varpath2 data.varname{eqn} '_' num2str(nstep) '.png'];
        
        figure(eqn)
        %fig = figure('OuterPosition',[50 50 900 900]);
        set(gcf,'Color','w')
        
    end

%end

DGsnapshot(:,:,:) = DGsolution(:,:,:);

for eqn = 1:data.Neqns 
    figure(data.fignum{eqn})
    clf 
    hold on
end
    
v1plot = zeros(data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));
v2plot = zeros(data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));
DGplot = zeros(data.Neqns,data.Nv1*length(data.plotlocs),data.Nv2*length(data.plotlocs));

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
        DGplot(:,Ixii+ixii,Ieta+ieta) = temp;
    end
    end    
end
end

for eqn = 1:data.Neqns
    figure(data.fignum{eqn})
    zplot(:,:) = DGplot(eqn,:,:);
    hold on
    surf(v1plot,v2plot,zplot)
    view(3)
    axis([data.boundsv1v2v3(1,1:4) data.boundsq(eqn,:) data.boundsq(eqn,:)])
    colorbar        
end
    

for eqn = 1:data.Neqns
    figure(data.fignum{eqn})
    title([data.varname{eqn} ' , t = ' num2str(tnow)],'FontSize',font)        
    shading interp
    xlabel(['x_1 N = ' num2str(data.Nv1)],'FontSize',font)
    ylabel(['x_2 N = ' num2str(data.Nv2)],'FontSize',font)
    colorbar
    set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %{
    if nstep == 1
        imwrite(imind,cm,data.gifname{eqn},'gif','Loopcount',inf,'DelayTime',delay);
    else%if ~mod(nstep,frameskip)
        imwrite(imind,cm,data.gifname{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    %}
    print(data.gifname{eqn},'-dpng')
       
    view(2)
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %{
    if nstep == 1
        imwrite(imind,cm,data.gifname2{eqn},'gif','Loopcount',inf,'DelayTime',delay);
    else%if ~mod(n,frameskip)
        imwrite(imind,cm,data.gifname2{eqn},'gif','WriteMode','append','DelayTime',delay);
    end
    %}
    print(data.gifname2{eqn},'-dpng')
    
end
end