%Experimental data
load measurements

%% Simulation for one model
[t, x] = Model0sim_ownpar();
Xsim = x(:,2);
Ssim = x(:,3);
Csim = x(:,4);

% Using curve fitter to interpolate the values
Xfit = createFit(t,Xsim);
Sfit  = createFit(t,Ssim);
Cfit = createFit(t,Csim);

%getting the values of the model (simulation) at the times of measurement
Ssim_int = Sfit(tS);
Xsim_int = Xfit(tX);
Csim_int = Cfit(tCO2);

%plotting
plot(tX,Xsim_int,"blue.--") %plotting simulation at measurement times
hold on
plot(tX,X,"blue.-")
plot(tS, Ssim_int, "red.--")
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

