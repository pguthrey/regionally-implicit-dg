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
problemname = 'ccadvection1D';
problemname = 'compressibleeuler1D_movingblob';
problemname = 'ccadvection1D_nonlinear';
problemname = 'burgers2D';



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
limiter = 'guthreylimiter';
limiter = 'rossmanithlimiter2';
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

Nvals = 20;%20 40 80]%[20 40 80]%;[40 80]
for M = 1
for r = 1
for basispick = 1%[1 2]
for methodi = 1%[1 2]
for iN = 1:length(Nvals)
    
clear data 
N = Nvals(iN);
%data.problemname = 'ccadvectiondiscont2D';
data.problemname = problemname;
[data] = problem_get_parameters(data);

data.maxiters = M;
data.residTol = 1e-06;
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

switch basispick
    case 1
        data.predictorbasis = 'Q';
    case 2
        data.predictorbasis = 'P';
end
data.correctorbasis = 'P';

data.r_param = r; 
data.verbose = false;
data.plotIC = false; 
data.plotwhilerunning = true;
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
        disp(['RIDG ' data.methodtype ' ' data.predictorbasis ])
    else
        disp(['LIDG ' data.methodtype ' ' data.predictorbasis ])   
    end
end
disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ])

switch data.r_param 
    case 2
        results_r2(iN,1) = N;
        results_r2(iN,2) = data.runtime;
        results_r2(iN,3) = data.Rel_L1err(1);
        results_r2(iN,4) = data.Rel_L2err(1);
        results_r2(iN,5) = data.Rel_Lierr(1);

        results_filtered_r2(iN,1) = N;
        results_filtered_r2(iN,2) = data.runtime;
        results_filtered_r2(iN,3) = data.Rel_L1err_filtered;
        results_filtered_r2(iN,4) = data.Rel_L2err_filtered;
        results_filtered_r2(iN,5) = data.Rel_Lierr_filtered;
    case 1
        results_r1(iN,1) = N;
        results_r1(iN,2) = data.runtime;
        results_r1(iN,3) = data.Rel_L1err(1);
        results_r1(iN,4) = data.Rel_L2err(1);
        results_r1(iN,5) = data.Rel_Lierr(1);

        results_filtered_r1(iN,1) = N;
        results_filtered_r1(iN,2) = data.runtime;
        results_filtered_r1(iN,3) = data.Rel_L1err_filtered;
        results_filtered_r1(iN,4) = data.Rel_L2err_filtered;
        results_filtered_r1(iN,5) = data.Rel_Lierr_filtered;
        
        results = results_r1;

    case 0
        results_li(iN,1) = N;
        results_li(iN,2) = data.runtime;
        results_li(iN,3) = data.Rel_L1err(1);
        results_li(iN,4) = data.Rel_L2err(1);
        results_li(iN,5) = data.Rel_Lierr(1);

        results_filtered_li(iN,1) = N;
        results_filtered_li(iN,2) = data.runtime;
        results_filtered_li(iN,3) = data.Rel_L1err_filtered;
        results_filtered_li(iN,4) = data.Rel_L2err_filtered;
        results_filtered_li(iN,5) = data.Rel_Lierr_filtered;
        
        results = results_li;
end


end


disp(' ')
disp(' ')
disp(' ')


end
end
end
end

%{
figure(99)
loglog(results(:,1),results(1,[3]).*(results(1,1)./results(:,1)).^data.M,'-ok')
hold on
loglog(results(:,1),results(:,[3:5]))



%}







