function [DGsolution,auxiliary,data] = problem_begin_timestep(DGsolution,auxiliary,data)

N(:,:,:,:) = DGsolution(1:data.Ls,:,:,:);

pressure_ideal = data.kbTi*N;

pressure_hartree = compute_hartree(DGsolution,auxiliary.pressure_hartree,data);
 
auxiliary.pressure_hartree = pressure_hartree;

[pressure_correlation,dpdn_excess_cell] = compute_correlation(DGsolution,data);

%Computes the total pressure as a modal expansion 
auxiliary.total_pressure = pressure_ideal ...
                           + pressure_hartree ...
                           + pressure_correlation; 
auxiliary.dpdn =  dpdn_excess_cell;
auxiliary.dpdn(1,:,:,:) = auxiliary.dpdn(1,:,:,:) + data.kbTi;

auxiliary.total_viscocity = compute_viscocity(DGsolution,data);


                       