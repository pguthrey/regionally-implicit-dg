function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = false;
data.boundsv1v2v3 = [-.5 .5 -.5 .5 0 0];

data.Tfinal = .2;

%NonCCAdvection
data.boundsq  = [-.05 3.5
                 -2 2
                 -2 2
                 -.05 16];
data.Neqns = 4;
data.varname{1} = 'Mass';
data.varname{2} = 'x momentum';
data.varname{3} = 'y momentum';
data.varname{4} = 'Energy';

data.Nplotvars = 9;
data.plotname{1} = 'Mass';
data.plotname{2} = 'x momentum';
data.plotname{3} = 'y momentum';
data.plotname{4} = 'Energy';
data.plotname{5} = 'Pressure';
data.plotname{6} = 'Velocity x';
data.plotname{7} = 'Velocity y';
data.plotname{8} = 'Specific Internal Energy';
data.plotname{9} = 'Speed';

data.boundsplot  = [-.1 3.8
                 -8.1 8.1 
                 -8.1 8.1 
                 -.1 16.1
                 -.1 12.1
                 -3.1 3.1
                 -3.1 3.1
                 -.1 16.1
                 -.1 3.1];

 data.Nlimitvars = 8;

data.appdata.gamma = 5/3; 

data.Fspeedmax = 5;
data.Gspeedmax = 5;
data.Hspeedmax = 0;

end

