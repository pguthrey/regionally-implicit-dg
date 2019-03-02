function [ q ] = DGeval(DGcoeffs,quadpoint,data)
% written by Pierson Guthrey

switch data.space_dims 
    case 1
        q = DGeval_1D(DGcoeffs,quadpoint,data);
    case 2
        q = DGeval_2D(DGcoeffs,quadpoint,data);
    case 3
        q = DGeval_3D(DGcoeffs,quadpoint,data);
end
