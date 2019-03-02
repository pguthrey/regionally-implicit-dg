function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = false;
data.boundsv1v2v3 = [data.v1_lb data.v1_ub 0 0 0 0];

%NonCCAdvection
data.boundsq  = [-.05 .1];

data.Neqns = 1;
data.varname{1} = 'Electron Density';

data.Nplotvars = 1;
data.plotname{1} = 'Electron Density';

data.boundsplot  = [-.1 .1];

data.Fspeedmax = 5;
data.Gspeedmax = 0;
data.Hspeedmax = 0;

end

