function [data] = problem_get_parameters(data)

data.space_dims = 1;

data.measureerror = false;
data.boundsv1v2v3 = [0 1 0 0 0 0];

data.Tfinal = .8;

%NonCCAdvection
data.boundsq  = [-.1 2
                 -2 2
                 -2 2
                 -.1 3];
data.Neqns = 3;
data.varname{1} = 'Mass';
data.varname{2} = 'x momentum ';
data.varname{3} = 'Energy';
%data.formulation = 'strong';

data.appdata.gamma = 5/3; 

data.Fspeedmax = 3;
data.Gspeedmax = 0;
data.Hspeedmax = 0;

data.Nplotvars = 6;
data.plotname{1} = 'Mass';
data.plotname{2} = 'Velocity';
data.plotname{3} = 'Energy';
data.plotname{4} = 'Pressure';
data.plotname{5} = 'Internal Energy';
data.plotname{6} = 'Momentum';

data.boundsplot  = [-.01 1.05
                 -.5 .5
                 -.01 1.6
                 -.01 1.05
                 -.01 2.5
                 -.5 .5 ];
             
data.linearcase = false;

end

