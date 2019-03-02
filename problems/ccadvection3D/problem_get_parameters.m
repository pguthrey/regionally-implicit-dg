function [data] = problem_get_parameters(data)

data.space_dims = 3;

data.measureerror = true;
data.boundsv1v2v3 = [-1 1 -1 1 -1 1];
data.Tfinal = 2;

%NonCCAdvection
data.boundsq  = [-.05 1.05];
data.Neqns = 1;
data.varname{1} = 'Mass';

data.Nplotvars = 1;
data.plotname{1} = 'Mass';

data.boundsplot  = [-.05 1.05];

data.appdata.nuf = 1;
data.appdata.nug = 1;
data.appdata.nuh = 1;
%%{
data.Fspeedmax = 1;
data.Gspeedmax = 1;
data.Hspeedmax = 1;
%}
data.linearcase = 'constantcoefficient';

end

