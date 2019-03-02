function [DGfinalsoln,data] = RIrDG_method_linear(DGsolution_old,data)
%Run the DG data.method for the given inputs
%   Detailed explanation goes here

deltatx = data.deltav1/data.Fspeedmax*data.cfl;
deltaty = data.deltav2/data.Gspeedmax*data.cfl;
deltatz = data.deltav3/data.Hspeedmax*data.cfl;

%%{
if strcmp(data.limiter,'epslimiter')
    disp('line 10 in RIrDG_method_nostore')
    deltaty = data.MAXDT ;    
end
%}
deltat = min([deltatx deltaty deltatz]);

Nsteps = max([10 ceil(data.Tfinal/deltat)]);%take at least 10 time steps
deltat = data.Tfinal/Nsteps;
data.nuv1 = deltat/data.deltav1;
data.nuv2 = deltat/data.deltav2;
data.nuv3 = deltat/data.deltav3;

[update,update_index,Kmax] = initialize_linear_update(data);

%maincell = update(:,:,5,2,2)
%tailcell = update(:,:,8,2,2)
%tailcellr1 = update(:,:,18,2,2)
%maincellr1 = update(:,:,13,2,2)
%save('testing','maincell','tailcell','maincellr1','tailcellr1')

switch data.r_param
    case 0
        %maincell_li = update(:,:,5,3,3)
        %tailcell_li = update(:,:,8,3,3)
        %save('testing','maincell_li','tailcell_li')
    case 1
        %maincell_r1 = update(:,:,13,3,3)
        %tailcell_r1 = update(:,:,18,3,3)
        %save('maincell_r1','tailcell_r1')
end

%{
tailcell_r1 =
  -0.333333333333333
maincell_r1 =
   0.416666666666667
%}
%plot_filtered_solution(DG_initialconditions,data)

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

%{
for iv1 = 1:length(data.v1centers)
for iv2 = 1:length(data.v2centers)
    X(iv1,iv2) = data.v1centers(iv1);
    Y(iv1,iv2) = data.v2centers(iv2);
    quadpoint = [X(iv1,iv2) Y(iv1,iv2)];
    U(iv1,iv2) = problem_F_FluxJacobian(0,quadpoint,data.appdata);
    V(iv1,iv2) = problem_G_FluxJacobian(0,quadpoint,data.appdata);
end
end
figure
quiver(X,Y,U,V)
keyboard
%}

%h1 = %waitbar(0,'Running Experiment...','OuterPosition', [100 400 300 75]);
tnow = 0;
nstep = 2;
%We timestep until we have reached the final time or we cannot timestep
experiment_dt = [];
experiment_nstep = [];
experiment_time = [];
data.Time = 0;
tic;
DGsolution_new = NaN(data.theta,data.Nv1,data.Nv2,data.Nv3);
for nstep = 1:Nsteps
    for iv1 = 1:data.Nv1
    for iv2 = 1:data.Nv2
    for iv3 = 1:data.Nv3       
        cell = DGsolution_old(:,iv1,iv2,iv3);
        for update_k = 1:Kmax
            if ~isnan(update_index(1,update_k,iv1,iv2,iv3))
                iv1_tilde = update_index(1,update_k,iv1,iv2,iv3);
                iv2_tilde = update_index(2,update_k,iv1,iv2,iv3);
                iv3_tilde = update_index(3,update_k,iv1,iv2,iv3);
                U = update(:,:,update_k,iv1,iv2,iv3);
                cell = cell - U*DGsolution_old(:,iv1_tilde,iv2_tilde,iv3_tilde);
            end
        end
        DGsolution_new(:,iv1,iv2,iv3) = cell;            
    end
    end
    end
    
    if 0%data.filter 
        scalemat = diag(sqrt(1:2:(2*data.M-1)));
        for ik = 1:length(DGsolution_new(1,:))
            DGsoln_tilde(:,ik) = scalemat*DGsolution_new(:,ik);    
        end
        solnfiltered = postProcess(DGsoln_tilde', data.deltav1, 10 );

        quadlocs = [
            -0.973906529
            -0.865063367
            -0.679409568
            -0.433395394
             -0.148874339
            0.148874339
            0.433395394
            0.679409568
            0.865063367
            0.973906529];
        quadwgts =[
            0.066671344
            0.149451349       
            0.219086363
            0.269266719        
            0.295524225        
            0.295524225
            0.269266719        
            0.219086363
            0.149451349
            0.066671344];
        for iv1 = 1:data.Nv1
        for ell = 1:data.theta
            temp = 0;
            for iq = 1:length(quadlocs)
                quadpoint = quadlocs(iq);
                weight = quadwgts(iq);
                phiT = testfunction_phi(quadpoint,ell,data);
                f = solnfiltered(iq+10*(iv1-1));
                temp = temp + phiT*f*weight/2;
            end                
            DGsolution_new(ell,iv1) = temp;
        end
        end
    end
    
    experiment_dt = [experiment_dt,deltat];
    experiment_nstep = [experiment_nstep,nstep-1];
    experiment_time = [experiment_time,tnow];
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
    DGsolution_old = DGsolution_new;
    %done = tnow/data.Tfinal;
    %waitbar(done, h1);
end
data.runtime = toc;
if data.verbose
    disp('Experiment complete')
end
%delete(h1)

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

filename = ['soln_' data.problemname '_r'  data.r_param '_M' num2str(data.M) '_N1' num2str(data.Nv1) '_N2' num2str(data.Nv2) '_N3' num2str(data.Nv3)  '_' data.limiter '_' data.basis];
if exist(['results\' filename ],'dir')
    rmdir(['results\' filename ],'s');
end


end


