function [data] = problem_get_parameters(data)

data.space_dims = 1;


data.measureerror = true;
data.boundsv1v2v3 = [-1 1 0 0 0 0];

data.Tfinal = 2;

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

data.linearcase = 'constantcoefficient';

end

