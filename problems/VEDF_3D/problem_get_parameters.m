function [data] = problem_get_parameters(data)

data.space_dims = 3;

data.boundsv1v2v3 = [-1 1 -1 1 -1 1];

data.Tfinal = .4;

%NonCCAdvection
data.boundsq  = [-.1 1.1];
data.Neqns = 7;
data.varname{1} = 'Density';
data.varname{2} = 'Velocity x';
data.varname{3} = 'Velocity y';
data.varname{4} = 'Velocity z';
data.varname{5} = 'Memory x';
data.varname{6} = 'Memory y';
data.varname{7} = 'Memory z';

data.Fspeedmax = 1.05;
data.Gspeedmax = 1.05;
data.Hspeedmax = 1.05;

data.appdata = 0;

data.Nplotvars = 4;
data.varname{1} = 'Density';
data.varname{2} = 'Velocity x';
data.varname{3} = 'Velocity y';
data.varname{4} = 'Velocity z';

data.boundsplot  = [-.1 1.1];
             
data.linearcase = 'nonlinear';

data.kbTi = 100;
data.kbTe = 1;
data.Z = 1;

data.measureerror = false;

end

