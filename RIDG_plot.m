clc
clear all
close all force
%drawnow
%distcomp.feature( 'LocalUseMpiexec', false );
dbstop if error
%dbclear all
warning('off','all')
%opengl('save', 'software')
enn = 1;

methods{1} = 'implicit';
methods{2} = 'picard';

%problemname = 'compressibleeuler1D_movingblob';
problemname = 'burgers2D';

limiter = 'none';

addpath(['problems/' problemname])
addpath(['limiters/' limiter])
addpath('predictor')
addpath('corrector')
addpath('compute')
addpath('testfunctions')
addpath('initialize')
addpath('plotting')
addpath('debug')
addpath('filter')
addpath('smartsolver')
addpath('RKDG')

minrun = 100000;
maxrun = 1e-04;
minerror = 1;
maxerror = 1e-07;
axes('NextPlot','add');
ploti = 1;

Nvals = [20]%[80 160 240 320 400]% 240 480]
%Nvals = 160

for M = 2
for corrbasispick = 1
for iN = 1:length(Nvals)
    
clear data 
N = Nvals(iN);
%data.problemname = 'ccadvectiondiscont2D';
data.problemname = problemname;
[data] = problem_get_parameters(data);

data.limiter = limiter;
           
data.Nv1 = N;
if data.space_dims >= 2
    data.Nv2 = N;
    if data.space_dims>= 3
        data.Nv3 = N;
    else
        data.Nv3 = 1;
    end
else
    data.Nv2 = 1;
    data.Nv3 = 1; 
end

data.r_param = 9; 
data.num_workers = 4; 
data.basiscombine = 'full';
data.basis = 'canonical';
data.extra = 'nonlinconvergence1d';
data.M = M;

data.predictorbasis = 'P';

switch corrbasispick
    case 1
        data.correctorbasis = 'Q';
    case 2
        data.correctorbasis = 'P';
end

%data.r_param = r; 
data.verbose = false;
data.plotIC = false; 
data.plotwhilerunning = false;
data.plotfinal = true;
data.makegif_conserved = false;
data.filter = false;%
data.smartsolver = true;

data.usewaitbars = true;
data.check_symmetry = false;
data.savefile = false;
data.email = false;
data.store = false;
data.check_conservation = false;
data.flags(1:data.Nv1) = false;

%data.methodtype = methods{methodi};

%%{
%delete(gcp('nocreate'))
%parpool('local',data.num_workers)
%poolobj = gcp;

%{
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('local',data.num_workers);
end
%}
 

end


end
end
%}

for M = 2%[4 6]
for r = [1]
for methodi = [1]
for predbasispick = 1%[1 2]
for corrbasispick = 2%[2 1]
for maxiters = M%[M-1:M+2]
for residtol = 1e-06%[1e-04 1e-06 1e-08]
for iN = 1:length(Nvals)
    
clear data 
N = Nvals(iN);
%data.problemname = 'ccadvectiondiscont2D';
data.problemname = problemname;
[data] = problem_get_parameters(data);

data.residTol = residtol;
data.maxiters = maxiters; 
data.limiter = limiter;
           
data.Nv1 = N;
if data.space_dims >= 2
    data.Nv2 = N;
    if data.space_dims>= 3
        data.Nv3 = N;
    else
        data.Nv3 = 1;
    end
else
    data.Nv2 = 1;
    data.Nv3 = 1; 
end

data.num_workers = 4; 
data.basiscombine = 'full';
data.basis = 'canonical';
data.extra = 'nonlinconvergence1d';
data.M = M;
%data.predictor_solver = 'Fsolve';
data.predictor_solver = 'NewtonIteration';

switch predbasispick
    case 1
        data.predictorbasis = 'Q';
    case 2
        data.predictorbasis = 'P';
end

switch corrbasispick
    case 1
        data.correctorbasis = 'Q';
    case 2
        data.correctorbasis = 'P';
end
data.r_param = r; 
data.verbose = false;
data.plotIC = false; 
data.plotwhilerunning = false;
data.plotfinal = true;
data.makegif_conserved = false;
data.filter = false;%
data.smartsolver = true;

data.usewaitbars = true;
data.check_symmetry = false;
data.savefile = false;
data.email = false;
data.store = false;
data.check_conservation = false;
data.flags(1:data.Nv1) = false;

data.methodtype = methods{methodi};

%%{
%delete(gcp('nocreate'))
%parpool('local',data.num_workers)
%poolobj = gcp;
%{
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('local',data.num_workers);
end
%}

[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
[data] = initialize_gauss_quadrature(data,data.M);
[data] = initialize_innerproduct_quadrature(data);
[data] = initialize_mesh(data);
[data] = initialize_precompute_testfunctions(data);
[data] = initialize_formIntegrals(data,data.Tfinal);        
[data] = limiter_initialize(data);

DG_initialconditions = projection_DGL2_Proj(@(point) problem_IC(point,data.appdata),data);

switch data.linearcase
    case 'constantcoefficient'
    [DGfinalsoln,data] = RIrDG_method_cc(DG_initialconditions,data);        
    case 'linear'
    [DGfinalsoln,data] = RIrDG_method_linear(DG_initialconditions,data);
    case 'nonlinear'
    [DGfinalsoln,data] = RIrDG_method_nostore(DG_initialconditions,data);    
end

if data.measureerror 
    if data.filter
        [ data ] = finalize_error_filtered(DGfinalsoln,data);
        
        disp(['L1 relative error = ' num2str(data.Rel_L1err)])
        disp(['L2 relative error = ' num2str(data.Rel_L2err)])
        disp(['Linf relative error = ' num2str(data.Rel_Lierr)])
        disp(['L1 relative error (filtered) = ' num2str(data.Rel_L1err_filtered)])
        disp(['L2 relative error (filtered) = ' num2str(data.Rel_L2err_filtered)])
        disp(['Linf relative error (filtered) = ' num2str(data.Rel_Lierr_filtered)])
    else
        [ data ] = finalize_error(DGfinalsoln,data);
        %{
        disp(['L1 relative error = ' num2str(data.Rel_L1err)])
        disp(['L2 relative error = ' num2str(data.Rel_L2err)])
        disp(['Linf relative error = ' num2str(data.Rel_Lierr)])
        disp(['runtime = ' num2str(data.runtime)])
        %}
    end
else
    data.Rel_L1err = NaN;
    data.Rel_L2err = NaN;
    data.Rel_Lierr = NaN;
    data.Rel_L1err_filtered = NaN;
    data.Rel_L2err_filtered = NaN;
    data.Rel_Lierr_filtered = NaN;
end

if iN == 1
    if r == 1
        runtitle = ['RIDG ' data.methodtype ' ' data.predictorbasis ' ' data.correctorbasis ' ' num2str(data.cfl) ' ' num2str(data.maxiters) ' ' num2str(data.residTol)];
    else
        runtitle = ['LIDG ' data.methodtype ' ' data.predictorbasis ' ' data.correctorbasis ' ' num2str(data.cfl) ' ' num2str(data.maxiters) ' ' num2str(data.residTol)]   ;
    end    
    disp(runtitle)
end    
    
quality = log(1/(data.Rel_L1err(1)*data.runtime));

disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ' ' num2str(quality)])

error_list(iN,1) = data.Rel_L1err(1);
runtime_list(iN,1) = data.runtime;

end



end
end
end
end
end
end
end

%{
close all
x = 0:.01:2*pi;
for k = 1:length(x)
   y(k) = problem_solution(data.Tfinal,x(k),data) ;
end
plot(x,y)
%}













