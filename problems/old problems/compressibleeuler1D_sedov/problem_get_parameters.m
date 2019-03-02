function [data] = problem_get_parameters(data)

data.space_dims = 1;

data.measureerror = false;
data.boundsv1v2v3 = [-.5 .5 0 0 0 0];

data.Tfinal = .2;

%NonCCAdvection
data.boundsq  = [-.05 3.5
                 -2 2
                 -.05 16];
data.Neqns = 3;
data.varname{1} = 'Mass';
data.varname{2} = 'x momentum ';
data.varname{3} = 'Energy';
%data.formulation = 'strong';

data.Nplotvars = 6;
data.plotname{1} = 'Mass';
data.plotname{2} = 'Velocity';
data.plotname{3} = 'Energy';
data.plotname{4} = 'Pressure';
data.plotname{5} = 'Internal Energy';
data.plotname{6} = 'Momentum';

data.boundsplot  = [-.1 3.8
                 -2.1 2.1
                 -.1 16.1
                 -.1 12.1
                 -.1 16
                 -8.1 8.1 ];

data.appdata.gamma = 5/3; 

data.Fspeedmax = 3;
data.Gspeedmax = 0;
data.Hspeedmax = 0;

end

