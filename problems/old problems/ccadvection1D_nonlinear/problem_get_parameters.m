function [data] = problem_get_parameters(data)

data.space_dims = 1;

data.measureerror = true;
data.boundsv1v2v3 = [-1 1 0 0 0 0];

data.Tfinal = 4;

%NonCCAdvection
data.boundsq  = [-1.1 1.1];
data.Neqns = 1;
data.varname{1} = 'Color';
%data.formulation = 'strong';

data.appdata.nuf = 1;
data.appdata.nug = 0;
data.appdata.nuh = 0;
data.Fspeedmax = 1;
data.Gspeedmax = 0;
data.Hspeedmax = 0;

data.Nplotvars = 1;
data.plotname{1} = 'Color';
data.boundsplot  = [-1.1 1.1];

data.linearcase = 'nonlinear';


data.problem_F = @(q,appdata) appdata.nuf*q ;
data.problem_JF = @(q,appdata) appdata.nuf;
data.problem_speedF_star = @(q,appdata) abs(appdata.nuf);
data.problem_dspeeddq_star = @(q,appdata) 0;

end

