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

quadrature{1} = 'quadrature';

problemname = 'VEDF_3D';

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

disp('warning: recent changes to basis')

Nvals = [5];
for r = [1]
for extra_quad = [0]
for quadmethod = [1]
for M = [3]
for maxiters = [1]
predbasispickvals = [1];%
corrbasispickvals = [1];
for bspick = [1]
for predbasispick = predbasispickvals(bspick)
for corrbasispick = corrbasispickvals(bspick)
for residtol = [1e-04]  
for iN = 1:length(Nvals)
for option = 1    
    
clear data    


data.badflag = 0;
data.predictor_quadrature = quadrature{quadmethod};

data.option = option;
%zdata.problemname = 'ccadvectiondiscont2D';
data.problemname = problemname;
[data] = problem_get_parameters(data);

data.residTol = residtol;
data.maxiters = maxiters; 
data.limiter = limiter;
%extra_quad = 0;%2; 


N = Nvals(iN);    
          
data.Nv1 = N;
data.Nv2 = 1;
data.Nv3 = 1; 
if data.space_dims >= 2
    data.Nv2 = N;
    if data.space_dims>= 3
        data.Nv3 = N;
    end
end

data.num_workers = 4; 
data.basiscombine = 'full';
data.basis = 'canonical';
data.extra = 'nonlinconvergence1d';
data.M = M;
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
data.verbose = true;
data.plotIC = false; 
data.plotwhilerunning = false;
data.plotfinal = false;
data.makegif_conserved = false;
data.filter = false;%
data.smartsolver = false;

data.waitbar = false;
data.check_symmetry = false;
data.savefile = false;
data.email = false;
data.store = false;
data.check_conservation = false;
data.flags(1:data.Nv1) = false;

data.methodtype = 'implicit';

data.chaching = 0 ;

[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
[data] = initialize_gauss_quadrature(data,M+extra_quad);
[data] = initialize_innerproduct_quadrature(data);
[data] = initialize_mesh(data);
[data] = initialize_precompute_testfunctions(data);
[data] = initialize_formIntegrals(data,data.Tfinal);        
[data] = limiter_initialize(data);
data.chaching = 0;

DG_initialconditions = projection_DGL2_Proj(@(point) problem_IC(point,data.appdata),data);

switch data.linearcase
    case 'constantcoefficient'
    [DGfinalsoln,data] = RIrDG_method_cc(DG_initialconditions,data);        
    case 'linear'
    [DGfinalsoln,data] = RIrDG_method_linear(DG_initialconditions,data);
    case 'nonlinear'
    [DGfinalsoln,data] = RIrDG_method_nonlinear(DG_initialconditions,data);    
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
        [ data ] = finalize_error_save(DGfinalsoln,data);
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
        runtitle = ['RIDG ' num2str(data.M) ' ' data.predictor_quadrature ' ' data.predictorbasis ' ' data.correctorbasis ' ' num2str(data.cfl) ' ' num2str(data.maxiters) ' ' num2str(data.residTol)];
    else
        runtitle = ['LIDG ' num2str(data.M) ' ' data.predictor_quadrature ' ' data.predictorbasis ' ' data.correctorbasis ' ' num2str(data.cfl) ' ' num2str(data.maxiters) ' ' num2str(data.residTol)];
    end    
    disp(runtitle)
end    
    
quality = log(1/(data.Rel_L1err(1)*data.runtime));

disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.numsteps) ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ' ' num2str(quality) ' ' num2str(data.badflag)])


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
end
end
end
