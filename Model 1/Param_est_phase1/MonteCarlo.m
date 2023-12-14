clear
fn = 'MonteCarlo_res_261123_1825.xlsx';
% initializeExcelMC(fn) %run this when creating a new file

for i=1:100
    [init_par_bounds, loBounds, upBounds] = bounds(); %defining bounds for initial parameters
    initparams  = initialize(init_par_bounds);
    [results, EXITFLAG] = param_est_V1(initparams, loBounds, upBounds);
    writetoExcel(fn, results, EXITFLAG, initparams)
end

%% Function for defining bounds for initial parameters and returning arrays with upper and lower bounds
function [init_par_bounds, loBounds, upBounds] = bounds()
    %define bounds
    mu_max_bounds = [0 1];
    K_s_bounds = [0 5];
    k_d_bounds = [0 1];
    Y_xs_bounds = [0 1];
    Y_xco2_bounds = [0 1];
    
    init_par_bounds = [mu_max_bounds ; K_s_bounds ; k_d_bounds ; Y_xs_bounds ; Y_xco2_bounds];
    
    loBounds = [mu_max_bounds(1) K_s_bounds(1) k_d_bounds(1) Y_xs_bounds(1) Y_xco2_bounds(1)];
    upBounds = [mu_max_bounds(2) K_s_bounds(2) k_d_bounds(2) Y_xs_bounds(2) Y_xco2_bounds(2)];
end

%% Select random initial values within boundaries
function initparams = initialize(init_par_bounds)
    initparams = zeros(1,5);
    for i=1:length(init_par_bounds)
        initparams(i) = rand(1)*(init_par_bounds(i,2)-init_par_bounds(i,1))+init_par_bounds(i,1); % random number from uniform distribution within boundaries
    end
end

%% Append results from parameter estimation to Excel file
function writetoExcel(fn,results, exitflag, initparams)
    % writematrix(initparams,fn, 'WriteMode','append')
    writematrix([results exitflag initparams], fn, 'WriteMode','append')
    writematrix([], fn, 'WriteMode','append')
end

function initializeExcelMC(fn)
    writematrix(["mu_max" "ks" "kd" "Yxs" "Yxco2" "Exitflag" "Initial values"], fn, "WriteMode", "append")
    writematrix([], fn, 'WriteMode','append')
end

