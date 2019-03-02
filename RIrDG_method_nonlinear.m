function [DGfinalsoln,data] = RIrDG_method_nonlinear(DGsolution_old,data)
%Run the DG data.method for the given inputs
%   Detailed explanation goes here
               
deltatx = data.deltav1/data.Fspeedmax*data.cfl;
deltaty = data.deltav2/data.Gspeedmax*data.cfl;
deltatz = data.deltav3/data.Hspeedmax*data.cfl;

%%{
%}
deltat = min([deltatx deltaty deltatz]);
data.nuv1 = deltat/data.deltav1;
data.nuv2 = deltat/data.deltav2;
data.nuv3 = deltat/data.deltav3;

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
nstep = 0;
%We timestep until we have reached the final time or we cannot timestep
experiment_dt = [];
experiment_nstep = [];
experiment_time = [];
data.Time = 0;
check = 1;
deltatnew = deltat;
stepscore  = 0;

%plot_filtered_solution(DGsolution_old,data,0)

start_cputime = cputime;
tic;

[DGsolution_old,auxiliary,data] = problem_initialize_auxiliary(DGsolution_old,data);

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
        if data.linearcase
            [data] = initialize_formIntegrals(data,timeleft);
        end
    end    
    if data.verbose
        disp(['Attempting dt = ' num2str(deltat) ' at t = ' num2str(tnow)])
    end
    data.deltat = deltat;
    data.tnow = tnow;
    
    [DGsolution_old,auxiliary,data] = problem_begin_timestep(DGsolution_old,auxiliary,data);
    
    %Take one timestep    
    switch data.methodtype
        case 'implicit'
            [ DGprediction,auxiliary,predmaxF,predmaxG,predmaxH] = predictor_main_implicit(DGsolution_old,auxiliary,data);
            predmaxF = 0;
        otherwise
            error(' Method type not yet implemented. ')          
    end
        
    
    [DGprediction,auxiliary] = problem_post_prediction(DGprediction,auxiliary);

    [DGcorrected,auxiliary,corrmaxF,corrmaxG,corrmaxH,data] = corrector_main_noghosts(DGprediction,DGsolution_old,auxiliary,data); 
    
    [DGcorrected,auxiliary] = problem_post_correction(DGcorrected,auxiliary);

    
    data.deltat = deltat;
    [DGsolution_new,data] = limiter_phi(DGcorrected,DGsolution_old,data);
    
    experiment_dt = [experiment_dt,deltat];
    experiment_nstep = [experiment_nstep,nstep-1];
    experiment_time = [experiment_time,tnow];
    %Determine what our dt should have been for this step
    maxspeeds = [predmaxF corrmaxF];
    maxspeedx = max([predmaxF corrmaxF]);
    maxspeedy = max([predmaxG corrmaxG]);
    maxspeedz = max([predmaxH corrmaxH]);
    
    %{
    maxspeedx = 1.05;
    if data.space_dims >= 2
        maxspeedy = 1.05;
    end
      %}

  
    dtspeedF = data.deltav1/maxspeedx*data.cfl;
    dtspeedG = data.deltav2/maxspeedy*data.cfl;
    dtspeedH = data.deltav3/maxspeedz*data.cfl;
        
    deltatnew = min([dtspeedF dtspeedG dtspeedH]);
    if deltat < 1e-9
        if data.verbose
            disp('small timestep, ending experiment')
        end
        tnow = data.Tfinal;                
    elseif deltat > deltatnew
        %Reject timestep
        if data.verbose
            disp('reject step, trying smaller dt')
        end
        deltat = .98*deltatnew;
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;    
        data.nuv3 = deltat/data.deltav3;    
        data.badflag = data.badflag+1;
        disp(' I rejected a step! I cannot afford to do that')
        %keyboard
        
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
        nstep = nstep + 1;
        data.Time(nstep) = tnow ; 
        if data.makegif_conserved
            [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
        elseif data.plotwhilerunning 
            problem_plotting_whilerunning(DGsolution_new,tnow,data)
        end
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
        deltat = .98*deltatnew;
        data.nuv1 = deltat/data.deltav1;
        data.nuv2 = deltat/data.deltav2;
        data.nuv3 = deltat/data.deltav3;
        if data.makegif_conserved
            [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
        elseif data.plotwhilerunning 
            problem_plotting_whilerunning(DGsolution_new,tnow,data)
        end
        nstep = nstep + 1;
        data.Time(nstep) = tnow ;
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

data.Nt = nstep;

DGfinalsoln = DGsolution_new;

if data.makegif_conserved
    [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
elseif data.plotfinal  
    problem_plotting_whilerunning(DGsolution_new,tnow,data)
end

filename = ['soln_' data.problemname '_r'  data.r_param '_M' num2str(data.M) '_N1' num2str(data.Nv1) '_N2' num2str(data.Nv2) '_N3' num2str(data.Nv3)  '_' data.limiter '_' data.basis];

if exist(['results\' filename ],'dir')
    rmdir(['results\' filename ],'s');
end

end


