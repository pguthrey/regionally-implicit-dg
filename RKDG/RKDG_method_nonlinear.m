function [DGfinalsoln,data] = RKDG_method_nonlinear(DGsolution_old,data)
%Run the DG data.method for the given inputs
%   Detailed explanation goes here
           
DGsolution_IC = DGsolution_old;

%%{
if strcmp(data.limiter,'epslimiter')
    disp('line 10 in RIrDG_method_nostore')
    deltaty = data.MAXDT ; 
end
%}

if data.check_conservation
    [data] = check_conservation(DGsolution_old,data);
end
%{
switch data.limiter
    case 'none'
        DGsolution_old = DG_initialconditions;
    case 'DGlimiter'
        DGsolution_old = limiter_DGlimiter(DG_initialconditions,DG_initialconditions,data);
    otherwise
        error('Invalid limiter selected');
end
%}  
if data.makegif_conserved
    [data] = plotting_makegif(DGsolution_old,0,data,0);
elseif data.plotIC
    problem_plotting_IC(DGsolution_old,data) 
end        

if data.waitbar
    h1 = waitbar(0,'Running Experiment...','OuterPosition', [100 400 300 75]);
end
tnow = 0;
nstep = 2;
%We timestep until we have reached the final time or we cannot timestep
experiment_dt = [];
experiment_nstep = [];
experiment_time = [];
data.Time = 0;
check = 1;
stepscore  = 0;


%Obtain the tableau]
switch data.M
    case 1
        RKa = [];
        RKb = [1];
        RKc = [];
        RKp = 1;%number of stages
    case 2
        RKa = [1/2];
        RKb = [0 1];
        RKc = [1/2];
        RKp = 2;%number of stages
    case 4
        RKa =[[1/2 0 0];
              [0 1/2 0];
              [0 0 1]];
        RKb = [1/6 1/3 1/3 1/6];
        RKc = [1/2 1/2 1];
        RKp = 4;%number of stages
    otherwise
        error('gimme dat tableau')
end 

%pick your cfl restriction
switch data.space_dims
    case 1
        switch data.M
            case 1
                data.cfl = 1.0;
            case 2
                data.cfl = 0.35;
            case 4
                data.cfl = 0.1;        
            otherwise
                error('gimme dat tableau')
        end 
    case 2
        switch data.M
            case 1
                error('gimme dat cfl')
            case 2
                data.cfl = 0.1;                
            case 4
                data.cfl = 0.05;   
            otherwise
                error('gimme dat cfl')
        end         
    otherwise
        error('gimme dat cfl')        
end


deltatx = data.deltav1/data.Fspeedmax*data.cfl;
deltaty = data.deltav2/data.Gspeedmax*data.cfl;
deltatz = data.deltav3/data.Hspeedmax*data.cfl;
deltat = min([deltatx deltaty deltatz]);
data.nuv1 = deltat/data.deltav1;
data.nuv2 = deltat/data.deltav2;
data.nuv3 = deltat/data.deltav3;
deltatnew = deltat;

%plot_filtered_solution(DGsolution_old,data,0)

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;


start_cputime = cputime;
tic;

RKh = 2;

v1centers = data.v1centers;
v2centers = data.v2centers;
v3centers = data.v3centers;

wrapind = @(x, n) (1 + mod(x-1, n));

global_speedmaxF = 0;
global_speedmaxG = 0;
global_speedmaxH = 0;

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;

local_maxspeedG = 0;
local_maxspeedH = 0;

space_dims = data.space_dims;

PHIinv = data.PHI_norm_inv;





while check 
    timeleft = data.Tfinal - tnow;
    if deltat >= timeleft
        %Reduce timestep
        if data.verbose
            disp('final step, trying time left as deltat')
        end
        deltat = timeleft;
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;
        data.nuv3 = deltat/data.deltav3;
    end    
    if data.verbose
        disp(['Attempting dt = ' num2str(deltat) ' at t = ' num2str(tnow)])
    end
    
    
    global_speedmaxF = 0;
    global_speedmaxG = 0;
    global_speedmaxH = 0;
          
    
    DGstagein = zeros(data.theta,Nv1,Nv2,Nv3);
    RK_k = zeros(data.theta,length(RKb),Nv1,Nv2,Nv3);
    
    for stage = 1:length(RKb)
        for iv1 = 1:Nv1
        for iv2 = 1:Nv2
        for iv3 = 1:Nv3        
            DGstagein(:,iv1,iv2,iv3) = DGsolution_old(:,iv1,iv2,iv3);
            for ind = 1:(stage-1)
                DGstagein(:,iv1,iv2,iv3) = DGstagein(:,iv1,iv2,iv3) + RKh*RK_k(:,ind,iv1,iv2,iv3)*RKa(stage-1,ind);
            end            
        end
        end
        end
            
        for iv1 = 1:Nv1
        for iv2 = 1:Nv2
        for iv3 = 1:Nv3
                v1center = v1centers(iv1); 
                v2center = v2centers(iv2); 
                v3center = v3centers(iv3); 
                cellcenter = [v1center v2center v3center];

                qstar = DGstagein(:,iv1,iv2,iv3);
                qeast = DGstagein(:,wrapind(iv1+1,Nv1),iv2,iv3);
                qwest = DGstagein(:,wrapind(iv1-1,Nv1),iv2,iv3);       
                [rhs,local_maxspeedF] = problem_RKDG_local_residual_eastwest(qstar,qeast,qwest,data,cellcenter);
                [rhs] = problem_RKDG_local_source(qstar,rhs,data,cellcenter);
                if space_dims >= 2
                    qnort = DGstagein(:,iv1,wrapind(iv2+1,Nv2),iv3);
                    qsout = DGstagein(:,iv1,wrapind(iv2-1,Nv2),iv3);
                    [rhs,local_maxspeedG] = problem_RKDG_local_residual_nortsout(qstar,qnort,qsout,rhs,data,cellcenter);
                    if space_dims >= 3
                        quppr = DGstagein(:,iv1,iv2,wrapind(iv3+1,Nv3));
                        qdown = DGstagein(:,iv1,iv2,wrapind(iv3-1,Nv3));
                        [rhs,local_maxspeedG] = problem_RKDG_local_residual_upprdown(qstar,quppr,qdown,rhs,data,cellcenter);
                    end
                end
                RK_k(:,stage,iv1,iv2,iv3) = -PHIinv*rhs;
                global_speedmaxF = max([local_maxspeedF global_speedmaxF]);
                global_speedmaxG = max([local_maxspeedG global_speedmaxG]);
                global_speedmaxH = max([local_maxspeedH global_speedmaxH]);       
        end
        end
        end
    end    
    
    maxspeedx = global_speedmaxF;
    maxspeedy = global_speedmaxG;
    maxspeedz = global_speedmaxH;
    
    DGsolution_new = DGsolution_old;        
    for iv1 = 1:Nv1
    for iv2 = 1:Nv2
    for iv3 = 1:Nv3
        for ind = 1:length(RKb)
            DGsolution_new(:,iv1,iv2,iv3) = DGsolution_new(:,iv1,iv2,iv3) + RKh*RK_k(:,ind,iv1,iv2,iv3)*RKb(ind);            
        end
    end
    end
    end
       %{
        

        
        
        
            [RK_k(:,1,iv1,iv2,iv3),speedmaxF,speedmaxG,speedmaxH] = RKDG_global_rhs(DGsolution_old,data);
            qstar = DGprevious(:,iv1,iv2,iv3);
            qeast = DGprevious(:,wrapind(iv1+1,Nv1),iv2,iv3);
            qwest = DGprevious(:,wrapind(iv1-1,Nv1),iv2,iv3);       
            [rhs,local_maxspeedF] = problem_RKDG_local_residual_eastwest(qstar,qeast,qwest,data,cellcenter);
            if space_dims >= 2
                qnort = DGprevious(:,iv1,wrapind(iv2+1,Nv2),iv3);
                qsout = DGprevious(:,iv1,wrapind(iv2-1,Nv2),iv3);
                [rhs,local_maxspeedG] = problem_RKDG_local_residual_nortsout(qstar,qnort,qsout,rhs,data,cellcenter);
            end
            global_speedmaxF = max([local_maxspeedF global_speedmaxF]);
            global_speedmaxG = max([local_maxspeedG global_speedmaxG]);
            global_speedmaxH = max([local_maxspeedH global_speedmaxH]);       
            rhs_global(:,:,iv1,iv2,iv3) = -PHIinv*rhs;
    end
    end
    end    
    
    
    

    switch data.space_dims
        case 1
            %space dimensions 1
            
            [RK_k(:,:,1),speedmaxF,speedmaxG,speedmaxH] = RKDG_global_rhs(DGsolution_old,data);
            maxspeedx = max(maxspeedx,speedmaxF);
            maxspeedy = max(maxspeedy,speedmaxG);
            maxspeedz = max(maxspeedz,speedmaxH);
            
            
            for stage = 2:length(RKb)
                DGstagein = DGsolution_old;
                for ind = 1:(stage-1)
                    DGstagein = DGstagein + RKh*RK_k(:,:,ind)*RKa(stage-1,ind);
                end
                [RK_k(:,:,stage),speedmaxF,speedmaxG,speedmaxH] = RKDG_global_rhs(DGstagein,data);
                maxspeedx = max(maxspeedx,speedmaxF);
                maxspeedy = max(maxspeedy,speedmaxG);
                maxspeedz = max(maxspeedz,speedmaxH);
            end

            DGsolution_new = DGsolution_old;
            for ind = 1:length(RKb)
                for iv1 = 1:Nv1
                    DGsolution_new(:,iv1) = DGsolution_new(:,iv1) + RKh*RK_k(:,iv1,ind)*RKb(ind);
                end
            end 
            
        case 2
            %space dimensions 2
            [RK_k(:,:,:,1),speedmaxF,speedmaxG,speedmaxH] = RKDG_global_rhs(DGsolution_old,data);
            maxspeedx = max(maxspeedx,speedmaxF);
            maxspeedy = max(maxspeedy,speedmaxG);
            maxspeedz = max(maxspeedz,speedmaxH);
            for stage = 2:length(RKb)
                DGstagein = DGsolution_old;
                for ind = 1:(stage-1)
                    DGstagein = DGstagein + RKh*RK_k(:,:,:,ind)*RKa(stage-1,ind);
                end
                [RK_k(:,:,:,stage),speedmaxF,speedmaxG,speedmaxH] = RKDG_global_rhs(DGstagein,data);
                maxspeedx = max(maxspeedx,speedmaxF);
                maxspeedy = max(maxspeedy,speedmaxG);
                maxspeedz = max(maxspeedz,speedmaxH);
            end

            DGsolution_new = DGsolution_old;
            for ind = 1:length(RKb)
                for iv1 = 1:Nv1
                    for iv2 = 1:Nv2
                        DGsolution_new(:,iv1,iv2) = DGsolution_new(:,iv1,iv2) + RKh*RK_k(:,iv1,iv2,ind)*RKb(ind);
                    end
                end
            end
            
    end
    %}

       
    %data.deltat = deltat;
    %[DGsolution_new,data] = limiter_phi(DGcorrected,DGsolution_old,data);
    
    
    
    experiment_dt = [experiment_dt,deltat];
    experiment_nstep = [experiment_nstep,nstep-1];
    experiment_time = [experiment_time,tnow];
    %Determine what our dt should have been for this step
    
    dtspeedF = data.deltav1/maxspeedx*data.cfl;
    dtspeedG = data.deltav2/maxspeedy*data.cfl;
    dtspeedH = data.deltav3/maxspeedz*data.cfl;
    
    deltatnew = min([dtspeedF dtspeedG dtspeedH]);
    if deltat > deltatnew
        %Reject timestep
        if data.verbose
            disp('reject step, trying smaller dt')
        end
        deltat = .95*deltatnew;% deltatnew*.95;
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;
        data.nuv3 = deltat/data.deltav3;
    elseif (deltat == deltatnew)
        %Accept Timestep, no changes to dt
        stepscore = 0;
        if data.verbose
            disp('accept, keep same dt')
        end
        if data.check_conservation
            [data] = check_conservation(DGsolution_new,data);
        end
        tnow = tnow + deltat;
        data.Time(nstep) = tnow ;
        if data.makegif_conserved
            [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
        elseif data.plotwhilerunning 
            problem_plotting_whilerunning(DGsolution_new,tnow,data)
        end
        nstep = nstep + 1;
        DGsolution_old = DGsolution_new;
        if data.waitbar
            done = tnow/data.Tfinal;
            waitbar(done, h1);
        end
    else 
        %Accept timestep, increase dt
        if data.verbose
            disp('accept, change dt')
        end
        if data.check_conservation
            [data] = check_conservation(DGsolution_new,data);
        end        
        tnow = tnow + deltat;
        deltat = .95*deltatnew;
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;
        data.nuv3 = deltat/data.deltav3;
        data.Time(nstep) = tnow ;
        if data.makegif_conserved
            [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
        elseif data.plotwhilerunning 
            problem_plotting_whilerunning(DGsolution_new,tnow,data)
        end
        nstep = nstep + 1;
        DGsolution_old = DGsolution_new;
        if data.waitbar
            done = tnow/data.Tfinal;
            waitbar(done, h1);
        end
    end
    if data.verbose
        disp('------------------------------------')
    end
    check = (tnow < data.Tfinal);
end
data.numsteps = nstep;
data.runtime = toc;
data.cpu_runtime = cputime - start_cputime;
if data.verbose
    disp('Experiment complete')
end

if data.waitbar
    delete(h1)
end

%{
font = 16;
foldername = [data.problemname '_' num2str(data.Nv1) '_M' num2str(data.M) '_r' num2str(data.r_param)];
mkdir('results/',foldername)
path = ['results/' foldername '/'];
varpath = 'timestepping';
imgname = [path varpath '.png'];
figure
set(gcf,'Color','w')
clf 
hold on
semilogy(experiment_nstep,experiment_dt,'-o')
title('Step size at each time step','FontSize',font)        
xlabel('Time step, n','FontSize',font)
ylabel('log timestep size, \Delta t','FontSize',font)
set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
print(imgname,'-dpng')
mkdir('results/',data.problemname)
path = ['results/' foldername '/'];
varpath = 'timestepping2';
imgname = [path varpath '.png'];
figure
set(gcf,'Color','w')
clf 
hold on
plot(experiment_nstep,experiment_time,'-o')
title('Time t at each time step','FontSize',font)        
xlabel('Time step, n','FontSize',font)
ylabel('Time, t','FontSize',font)
%axis([0 max(experiment_nstep) 0 data.Tfinal])
set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
print(imgname,'-dpng')
mkdir('results/',data.problemname)
path = ['results/' foldername '/'];
varpath = 'timestepping3';
imgname = [path varpath '.png'];
figure
set(gcf,'Color','w')
clf 
hold on
semilogy(experiment_time,experiment_dt,'-o')
title('Step size vs time t','FontSize',font)        
%axis([0 data.Tfinal min(experiment_dt) max(experiment_dt)])
xlabel('Time, t','FontSize',font)
ylabel('log timestep size, \Delta t','FontSize',font)
set(gca,'FontSize',font)%,'XTick',-1:.25:1,'YTick',-1:.25:1)
print(imgname,'-dpng')    
%}

data.Nt = nstep-1;

DGfinalsoln = DGsolution_new;

if data.makegif_conserved
    [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
elseif data.plotfinal  
    problem_plotting_whilerunning(DGsolution_new,tnow,data)
end

%{
filename = ['soln_' data.problemname '_r'  data.r_param '_M' num2str(data.M) '_N1' num2str(data.Nv1) '_N2' num2str(data.Nv2) '_N3' num2str(data.Nv3)  '_' data.limiter '_' data.basis];

if exist(['results\' filename ],'dir')
    rmdir(['results\' filename ],'s');
end

%mkdir(['results\' filename ]);
%path = ['results\' filename '\'];
%solnfile =[path filename '.mat'];
%save(solnfile,'DGsolution','data')
%plot_filtered_solution(DGfinalsoln,data,data.Tfinal)

%}

end


