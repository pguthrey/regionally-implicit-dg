function [data] = plotting_makegif(DGsolution,tnow,data,nstep)
% written by Pierson Guthrey

switch data.space_dims
    case 1
        data = plotting_makegif_1D(DGsolution,tnow,data,nstep);
    case 2
        data = plotting_makegif_2D(DGsolution,tnow,data,nstep);
        %data = plotting_makegif_2Dto1D(DGsolution,tnow,data,nstep);
end


