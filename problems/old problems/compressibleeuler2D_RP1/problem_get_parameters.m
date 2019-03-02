function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = false;
data.boundsv1v2v3 = [0 1 0 1 0 0];

data.Tfinal = .8;

%NonCCAdvection
data.boundsq  = [-.1 1.8
                 -.8 .8
                 -.8 .8
                 -.05 3];
data.Neqns = 4;
data.varname{1} = 'Mass';
data.varname{2} = 'x momentum';
data.varname{3} = 'y momentum';
data.varname{4} = 'Energy';

data.Nplotvars = 9;
data.plotname{1} = 'Density';
data.plotname{2} = 'x momentum';
data.plotname{3} = 'y momentum';
data.plotname{4} = 'Energy';
data.plotname{5} = 'Pressure';
data.plotname{6} = 'Velocity x';
data.plotname{7} = 'Velocity y';
data.plotname{8} = 'Specific Internal Energy';
data.plotname{9} = 'Speed';

data.boundsplot  = [-.1 2
                 -.8 .8 
                 -.8 .8 
                 -.1 3.1
                 -.1 1.8
                 -1.5 1.7
                 -1.5 1.7
                 -.1 5.1
                 -0.1 3.1];

data.Nlimitvars = 8;

data.appdata.gamma = 5/3; 

data.Fspeedmax = 5;
data.Gspeedmax = 5;
data.Hspeedmax = 0;

end

