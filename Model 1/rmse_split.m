%% Simulation for split model
clear
%run simulation
[t, x] = Model1sim();
Xsim = x(:,2);
Ssim = x(:,3);
Csim = x(:,4);

%missing simulation between 18.3667 and 26.5 h
lb=21.3667;
ub=26.5;

%Experimental data
load measurements

%removing elements within the missing simulation (want to find RMSE only
%for the times where there are both measurements and simulation)
[tS, S] = removedIndices(tS,S,lb,ub);
[tX, X] = removedIndices(tX,X,lb,ub);
[tCO2, CO2m] = removedIndices(tCO2,CO2m,lb,ub);

% Using curve fitter to interpolate the values
Xfit = createFit(t,Xsim);
Sfit  = createFit(t,Ssim);
Cfit = createFit(t,Csim);

%getting the values of the model (simulation) at the times of measurement
Ssim_int = Sfit(tS);
Xsim_int = Xfit(tX);
Csim_int = Cfit(tCO2);

%plotting
figure()
plot(tX,Xsim_int,"blue.")
hold on
plot(t,Xsim,"blue--")
plot(tX,X,"blue.-")

plot(tS, Ssim_int, "red.")
plot(t,Ssim,"red--")
plot(tS,S,"red.-")

%calculating RMSE
RMSE_S = rmse(Ssim_int,S')
RMSE_X = rmse(Xsim_int,X')
RMSE_C = rmse(Csim_int,CO2m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult, gof] = createFit(t, state)
%Automatically generated by cftool
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, state );

% Set up fittype and options.
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'linear';
opts.Normalize = 'on';

% Fit sim to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'state vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 't', 'Interpreter', 'none' );
% ylabel( 'state', 'Interpreter', 'none' );
% grid on
end

function [time_state, state] = removedIndices(time, state, lb, ub)
    time_state = time;
    %removing elements within lb and ub
    removeindices = find(time >lb & time <ub);
    time_state(removeindices) = [];
    state(removeindices) = [];
end