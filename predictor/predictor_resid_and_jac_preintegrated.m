function [Jac_east_cell_temp,Jac_west_cell_temp, Jac_east_other_temp,Jac_west_other_temp, Jac_trunc_east_temp,Jac_trunc_west_temp,Jac_cell,trunceast_temp,truncwest_temp,fluxeast_temp,fluxwest_temp,residual_cell,maxspeedF] = predictor_resid_and_jac_preintegrated(qstar,qeast,qwest,qpast,data,cellcenter)
% written by Pierson Guthrey

if strcmp(data.predictorbasis,'Q')
    switch data.M
        case 4
            [Jac_east_cell_temp,Jac_west_cell_temp,...
                Jac_east_other_temp,Jac_west_other_temp, ...
                Jac_trunc_east_temp,Jac_trunc_west_temp, ...
                Jac_cell,trunceast_temp,...
                truncwest_temp,fluxeast_temp,fluxwest_temp,...
                residual_cell,maxspeedF] = problem_resid_and_jac_preintegrated_M4_Q(qstar,qeast,qwest,qpast,data,cellcenter);
            
            %[Jac_east_cell_temp,Jac_west_cell_temp, Jac_east_other_temp,Jac_west_other_temp, Jac_trunc_east_temp,Jac_trunc_west_temp,Jac_cell,trunceast_temp,truncwest_temp,fluxeast_temp,fluxwest_temp,residual_cell,maxspeedF] = problem_resid_and_jac_preintegrated_M4_Q(qstar,qeast,qwest,qpast,data,cellcenter);
        case 2
            [Jac_east_cell_temp,Jac_west_cell_temp, Jac_east_other_temp,Jac_west_other_temp, Jac_trunc_east_temp,Jac_trunc_west_temp,Jac_cell,trunceast_temp,truncwest_temp,fluxeast_temp,fluxwest_temp,residual_cell,maxspeedF] = problem_resid_and_jac_preintegrated_M2_Q(qstar,qeast,qwest,qpast,data,cellcenter);
        otherwise
            error('whoops')                    
    end
else
    switch data.M
        case 4
            [Jac_east_cell_temp,Jac_west_cell_temp,...
                Jac_east_other_temp,Jac_west_other_temp, ...
                Jac_trunc_east_temp,Jac_trunc_west_temp, ...
                Jac_cell,trunceast_temp,...
                truncwest_temp,fluxeast_temp,fluxwest_temp,...
                residual_cell,maxspeedF] = problem_resid_and_jac_preintegrated_M4_P(qstar,qeast,qwest,qpast,data,cellcenter);            
        otherwise
            error('whoops')                    
    end
    
end

