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
methods{2} = 'explicit';

%problemname = 'compressibleeuler2D_RP1';
%problemname = 'compressibleeuler2D_1direction';
%problemname = 'ccadvectiondiscont2D';
problemname = 'compressibleeuler1D_shocktube';
problemname = 'ccadvection1D_testlimiter';
problemname = 'ccadvection1D_dumbser_limiter';
problemname = 'ccadvection1D_discontinuous';
problemname = 'compressibleeuler1D_sedov';
problemname = 'ccadvection1D_discontinuous';
problemname = 'compressibleeuler2D_sedov';
problemname = 'compressibleeuler2D_RP1';
problemname = 'compressibleeuler2D_sedov_extreme';
problemname = 'ccadvection1D_discontinuous';
problemname = 'compressibleeuler2D_sedov';
problemname = 'compressibleeuler1D_sedov';
problemname = 'compressibleeuler1D_shocktube';
problemname = 'ccadvection1D_discontinuous';
problemname = 'compressibleeuler1D_shocktube';
problemname = 'ccadvection1D_alltest';
problemname = 'compressibleeuler1D_shocktube';
problemname = 'nonccadvection2D';
problemname = 'ccadvection1D';
problemname = 'nonccadvection2D';
problemname = 'ccadvection2D';
problemname = 'ccadvection1D_nonlinear';
problemname = 'compressibleeuler1D_movingblob';

limiter = 'guthrey_highmoments';
limiter = 'dumbserlimiter2';
limiter = 'dumbserlimiter';
limiter = 'dgweno';
limiter = 'rossmanithlimiter2';
limiter = 'subcell_WENO';
limiter = 'rossmanithlimiter';
limiter = 'rossmanithlimiterrepeater';
limiter = 'dumbserlimiter';
limiter = 'subcell_WENO';
limiter = 'subcell_WENO_region';
limiter = 'rossmanithlimiter2';
limiter = 'vertexlimiter';
limiter = 'subcell_aderweno';
limiter = 'rossmanithlimiter2';
limiter = 'guthreylimiter';
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

Nvals = [20 40 80]
for iN = 1:length(Nvals)
for M = 4
for r = 1
clear data 
N = Nvals(iN);
%data.problemname = 'ccadvectiondiscont2D';
data.problemname = problemname;
[data] = problem_get_parameters(data);

data.residTol = [10*eps];
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
data.predictorbasis = 'Q';
data.correctorbasis = 'P';

data.r_param = r; 
data.verbose = true;
data.plotIC = false; 
data.plotwhilerunning = false;%true;
data.plotfinal = true;
data.makegif_conserved = false;%true;
data.filter = false;%
data.smartsolver = true;

data.usewaitbars = true;
data.check_symmetry = false;
data.savefile = false;
data.email = false;
data.store = false;
data.check_conservation = false;
data.flags(1:data.Nv1) = false;

for methodi = 1
data.methodtype = methods{methodi};

close all force
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
[DGfinalsoln,data] = RKDG_method_nostore(DG_initialconditions,data);        

if data.measureerror 
    if data.filter
        [ data ] = finalize_error_filtered(DGfinalsoln,data);
        disp(['L1 relative error = ' num2str(data.Rel_L1err)])
        disp(['L2 relative error = ' num2str(data.Rel_L2err)])
        disp(['Linf relative error = ' num2str(data.Rel_Lierr)])
        disp(['L1 relative error (filtered) = ' num2str(data.Rel_L1err_filtered)])
        disp(['L2 relative error (filtered) = ' num2str(data.Rel_L2err_filtered)])
        disp(['Linf relative error (filtered) = ' num2str(data.Rel_Lierr_filtered)])
        disp(['Runtime = ' num2str(data.runtime)])
    else
        [ data ] = finalize_error(DGfinalsoln,data);
        disp(['L1 relative error = ' num2str(data.Rel_L1err)])
        disp(['L2 relative error = ' num2str(data.Rel_L2err)])
        disp(['Linf relative error = ' num2str(data.Rel_Lierr)])
        disp(['Runtime = ' num2str(data.runtime)])
    end
else
    data.Rel_L1err = NaN;
    data.Rel_L2err = NaN;
    data.Rel_Lierr = NaN;
    data.Rel_L1err_filtered = NaN;
    data.Rel_L2err_filtered = NaN;
    data.Rel_Lierr_filtered = NaN;
end

x = linspace(data.boundsv1v2v3(1),data.boundsv1v2v3(2),1000);

%{
disp('-------------------------')
disp(data.Nv1)
disp(data.r_param)
disp(data.runtime)
%}


results(iN,1) = N;
results(iN,2) = data.runtime;
results(iN,3) = data.Rel_L1err(1);
results(iN,4) = data.Rel_L2err(1);
results(iN,5) = data.Rel_Lierr(1);

results_filtered_li(iN,1) = N;
results_filtered_li(iN,2) = data.runtime;
results_filtered_li(iN,3) = data.Rel_L1err_filtered;
results_filtered_li(iN,4) = data.Rel_L2err_filtered;
results_filtered_li(iN,5) = data.Rel_Lierr_filtered;



end
end
end
end

figure(99)
loglog(results(:,1),results(1,[3]).*(results(1,1)./results(:,1)).^data.M,'--k')
hold on
loglog(results(:,1),results(:,[3:5]))
