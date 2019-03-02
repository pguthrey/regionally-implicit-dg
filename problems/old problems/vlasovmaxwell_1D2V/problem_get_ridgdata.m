function [data] = problem_get_ridgdata(data)


data.limiter = 'rossmanithlimiter';

[data] = problem_get_parameters(data);
data.space_dims = 2;
data.measureerror = false;
data.boundsv1v2v3 = [data.v1_lb data.v1_ub data.v2_lb data.v2_ub 0 0];
data.Fspeedmax = 5;
data.Gspeedmax = 5;
data.Hspeedmax = 0;
data.Neqns = 1;

data.Nv1 = 101;
data.Nv2 = 21%3;
data.Nv3 = 1;

data.basis = 'canonical';
data.extra = 'allM';
data.predictor_solver = 'NewtonIteration';

data.r_param = 1; 
data.verbose = false;
data.plotIC = false; 
data.plotwhilerunning = false;
data.plotfinal = false;
data.makegif_conserved = false;
data.usewaitbars = true;
data.savefile = false;
data.email = false;
data.store = false;
data.check_conservation = false;

data.methodtype = 'implicit';

[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
[data] = initialize_gauss_quadrature(data);
[data] = initialize_innerproduct_quadrature(data);
[data] = initialize_mesh(data);
[data] = initialize_precompute_testfunctions(data);
[data] = initialize_formIntegrals(data);

end

