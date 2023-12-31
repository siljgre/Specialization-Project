clear
fn = 'MonteCarlo_500.xlsx';
% initializeExcelMC(fn) %run this when creating a new file

for i=1:500
    [init_par_bounds, loBounds, upBounds] = bounds(); %defining bounds for initial parameters
    initparams  = initialize(init_par_bounds);
    [results, EXITFLAG] = param_est(initparams, loBounds, upBounds);
    writetoExcel(fn, results, EXITFLAG, initparams)
end

%% Function for defining bounds for initial parameters and returning arrays with upper and lower bounds

% Estimating:
function [init_par_bounds, loBounds, upBounds] = bounds()
    init_par_bounds = [[5 50]; [0 1];[0 2]; [0 2];[0 2]; [1 10]; [10 50]];
    loBounds = init_par_bounds(:, 1);
    upBounds = init_par_bounds(:, 2);
end

%% Select random initial values within boundaries
function initparams = initialize(init_par_bounds)
    initparams = zeros(1,length(init_par_bounds));
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
    writematrix(["gamma" "c" "mu_maxs" "mu_maxx" "mu_maxc" "K_m1" "K_m2" "Exitflag" "Initial values"], fn, "WriteMode", "append")
    writematrix([], fn, 'WriteMode','append')
end

