function [update,upindexlist,upmeshlocslist] = initialize_constantcoefficient_update(data)
% written by Pierson Guthrey

switch data.space_dims
    case 1
        [update,upindexlist,upmeshlocslist,LHS2,predictor_update2,corrector_update2] = initialize_constantcoefficient_update_1D(data);
    case 2
        [update,upindexlist,upmeshlocslist,LHS2,predictor_update2,corrector_update2,mats] = initialize_constantcoefficient_update_2D(data);
    case 3
        [update,upindexlist,upmeshlocslist,LHS2,predictor_update2,corrector_update2] = initialize_constantcoefficient_update_3D(data);
end

        