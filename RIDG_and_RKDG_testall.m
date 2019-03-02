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
quadrature{2} = 'quadrature_free_volumes';
quadrature{3} = 'quadrature_free_all';

problemname = 'burgers1D_sourceterm';
problemname = 'burgers2D_quadraturefree_2';
problemname = 'burgers2D_sourceterm';
%problemname = 'burgers1D_quadfree';

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

deltar = -2:2;
%%{
Nvals = [39 52 65 77 91 105 151];

Nvals = [11 22 33 44 55 66 122];
for M = 4
for corrbasispick = [1]
for pick = 1   
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
 
data.cfl = NaN;
data.methodtype = 'null';
%[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
[data] = initialize_gauss_quadrature(data,data.M);
[data] = initialize_innerproduct_quadrature(data);
[data] = initialize_mesh(data);
[data] = initialize_precompute_testfunctions(data);
[data] = initialize_formIntegrals(data,data.Tfinal);        
[data] = limiter_initialize(data);

DG_initialconditions = projection_DGL2_Proj(@(point) problem_IC(point,data.appdata),data);

switch pick
    case 1
        [DGfinalsoln,data] = RKDG_method_nonlinear(DG_initialconditions,data);  
    case 2
        [DGfinalsoln,data] = RKDG_method_nostore(DG_initialconditions,data);  
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
        %{
        disp('error')
        [ data ] = finalize_error(DGfinalsoln,data);
%        %{
        disp(['L1 relative error = ' num2str(data.Rel_L1err)])
        disp(['L2 relative error = ' num2str(data.Rel_L2err)])
        disp(['Linf relative error = ' num2str(data.Rel_Lierr)])
        disp(['runtime = ' num2str(data.runtime)])
        %}
        %disp(' error installer ')
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

quality = log(1/(data.Rel_L1err(1)*data.runtime));

if iN == 1
    runtitle = ['RKDG ' data.M ' ' data.correctorbasis ' ' num2str(data.cfl) ];
    disp(runtitle)
    disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.numsteps)  ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ' ' num2str(quality)])
else
    disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.numsteps)  ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ' ' num2str(quality)])
end

error_list(iN,1) = data.Rel_L1err(1);
runtime_list(iN,1) = data.runtime;

end

disp(' ')
disp(' ')
disp(' ')

%${
plotHandles(ploti) = loglog( runtime_list, error_list ,'-o','DisplayName', runtitle,'LineWidth',3);
plotLabels{ploti} = runtitle;
txt2 = ['\leftarrow ' runtitle];
text(runtime_list(1), error_list(1),txt2)
legend(plotHandles, plotLabels);
ploti = ploti + 1;
%loglog( runtime_list, error_list ,'DisplayName', runtitle);
minrun = min([minrun ; runtime_list]);
maxrun = max([maxrun ; runtime_list]);
minerror = min([minerror ; error_list]);
maxerror = max([maxerror ; error_list]);

%axis([minrun/1.2 1.1*maxrun minerror/1.2 maxerror*1.1])
legend('show')
set(gca,'XScale','log','YScale','log')
drawnow

end
end
end
%}


%=============================================================================================
% RIDG Test section
%=============================================================================================
%{

Nvals = [67+deltar 78+deltar 89+deltar 100+deltar 112+deltar 123+deltar 134+deltar 145+deltar 156+deltar];
Nvals = [200 30:70];
for r = [1]
for extra_quad = [0]
for quadmethod = [3]
for M = [4 6 8]
for maxiters = [3]
predbasispickvals = [1];%
corrbasispickvals = [1];
for bspick = [1]
for predbasispick = predbasispickvals(bspick)
for corrbasispick = corrbasispickvals(bspick)
for residtol = [1e-04] %[1e-04 1e-06 1e-08 1e-010]
%{
if predbasispick == 1
    Nvals = [100 29 41 57 71 85 101 130 145 161];
else    
    Nvals = [100 33 42 51 60 68 77 86 95 104];
end
%}
%Nvals = 5:20

deltar = -2:2;
switch M
    case 4
        %Nvals = [11+deltar 22+deltar 33+deltar 44+deltar 56+deltar 67+deltar];
        Nvals = [100 39 52 65 77 91 105 119];        
    case 6
        %Nvals = [11+deltar 22+deltar 33+deltar 44+deltar];
        Nvals = [100 13 26 39 53 66];
    case 8
        %Nvals = [11+deltar 22+deltar 33+deltar 44+deltar];
        Nvals = [10 3 8];
        residtol = [1e-06];
        quadmethod = [1];
end
    
%{
if (data.space_dims == 1) && (data.M == 4)
    Nvals = [20 13  25 40 52 67 80 93];
end
if (data.space_dims == 1) && (data.M == 6)
    Nvals = [20 13 26 40 53 67 80 94];
end
if (data.space_dims == 1) && (data.M == 8)
    Nvals = [20 3 8];
end
%}  


for iN = 1:length(Nvals)
for option = 1    
    
clear data    

data.quadfree_jac_init = false;

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
        runtitle = ['LIDG ' num2str(data.M) ' ' data.predictor_quadrature ' ' data.predictorbasis ' ' data.correctorbasis ' ' num2str(data.cfl) ' ' num2str(data.maxiters) ' ' num2str(data.residTol)]   ;
    end    
    disp(runtitle)
end    
    
quality = log(1/(data.Rel_L1err(1)*data.runtime));

disp([num2str(N) ' ' num2str(data.runtime) ' ' num2str(data.numsteps) ' ' num2str(data.Rel_L1err(1)) ' ' num2str(data.Rel_L2err(1)) ' ' num2str(data.Rel_Lierr(1)) ' ' num2str(quality) ' ' num2str(data.badflag)])


error_list(iN,1) = data.Rel_L1err(1);
runtime_list(iN,1) = data.runtime;

end
end

disp(' ')
disp(' ')
disp(' ')

figure(30)
plotHandles(ploti) = loglog( runtime_list, error_list ,'-o','DisplayName', runtitle);
plotLabels{ploti} = runtitle;
txt2 = ['\leftarrow ' runtitle];
text(runtime_list(1), error_list(1),txt2)
%legend(plotHandles, plotLabels);
ploti = ploti + 1;
%loglog( runtime_list, error_list ,'DisplayName', runtitle);
minrun = min([minrun ; runtime_list]);
maxrun = max([maxrun ; runtime_list]);
minerror = min([minerror ; error_list]);
maxerror = max([maxerror ; error_list]);
%axis([minrun/1.2 1.1*maxrun minerror/1.2 maxerror*1.1])
%axis([10 7000 1e-04 1])
legend('show')
set(gca,'XScale','log','YScale','log')
drawnow

end
end
end
end
end
end
end
end
end
%}
%{
close all
x = 0:.01:2*pi;
for k = 1:length(x)
   y(k) = problem_solution(data.Tfinal,x(k),data) ;
end
plot(x,y)
%}
%orders = diff(log(error_list))./-diff(log(Nvals'))












