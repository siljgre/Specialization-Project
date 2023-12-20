% two sets of parameters

function [t, x] = ode45model()

load measurements
%using only the timepoints of measurement for substrate and biomass for CO2 since the CO2 measurements are continous
measured = load("XSCt.mat");
CO2m = measured.sol(:,3);
tCO2 = measured.tspan;

% Initialize initial conditions
x0 = init_cond();
% Initialize parameters
par = param1();
% Define timespan differential equations are simulated for
% tspan = [0 21.3667];
tspan = [0 21.3667];
% Initialize options (ALWAYS EMPTY!)
options = [];

% Run numerical solver
[t1,x1] = ode45(@diff_eq,tspan,x0,options,par);

% skipping the addition - parameters not known

%% Addition of substrate
% x0 = x1(end,:);
% par = param2();
% par.F_in = 0.5/length(tspan);
% tspan = [21.3667 26.5];
% [t2,x2] = ode45(@diff_eq,tspan,x0,options,par);

%% After addition
[Xexp, Sexp, Cexp] = expdataphase3();
x0 = [2; Xexp(1); Sexp(1); Cexp(1)];

par = param3(); %new parameters
tspan = [26.500 50];
[t3,x3] = ode45(@diff_eq,tspan,x0,options,par);

% Concatenation of arrays
x = [x1; x3];
t = [t1; t3];

Biomass = x(:, 2);
Substrate = x(:, 3);
CO2 = x(:, 4);

% Plot figure
[tX, tS, CO2m, tCO2, S, X] = expdata();

figure(1)
plot(t, Biomass,"Color","blue")
hold on
plot(t, Substrate,"Color","red")
plot(t, CO2, "Color", "#77AC30")
plot(tX,X,"o","Color","blue")
plot(tS,S,"o","Color","red")
plot(tCO2,CO2m,"o","Color", "#77AC30")
rectangle('Position', [t1(end), 0 ,t3(1) - t1(end) , 25],'FaceColor', '#F0F0F0', 'LineStyle','--' );
axis([0 50 0 25]);
hold off
xlabel("Time [h]")
% legend("Biomass [g/L]", "Glucose [g/L]", "CO2 [%]")
legend("X [g/L]", "S [g/L]", "CO_2 [%]", "Measured X [g/L]", "Measured S [g/L]", "Measured CO_2 [%]")
end

% Function to define your differential equation
function dxdt = diff_eq(t,x,par)
Volume = x(1);
Biomass = x(2);
Substrate = x(3);
CO2 = x(4);

Volume(Volume<0.00001) = 0.00001; 
Biomass(Biomass<0) = 0;
Substrate(Substrate<0) = 0;
CO2(CO2<0) = 0;

% Stop emptying the tank when the tank is empty, to ensure that volume does
% not reach negative values
if Volume == 0.00001
    par.F_out = 0;
end

% Define differential equation
Volume_dot = par.F_in - par.F_out;
Biomass_dot = -(par.F_in./Volume).*Biomass+par.mu_max.*(Substrate./(par.K_s+Substrate)).*Biomass-par.k_d.*Biomass;
Substrate_dot = (par.F_in./Volume).*(par.S_in-Substrate)-(par.mu_max).*Substrate./((par.K_s+Substrate)).*(Biomass./par.Y_xs);
CO2_dot = par.mu_max.*(Substrate./(par.K_s+Substrate)).*(Biomass./par.Y_xco2)-par.q_air.*CO2;

% Save differential equation in output vector
dxdt = [Volume_dot; Biomass_dot; Substrate_dot; CO2_dot];
end

% Function to define all your constant parameters
function par = param1()
%%%%%%%%%%%%%%%
par.S_in = 100;
par.F_in = 0;
par.F_out = 0;
par.q_air = 2;
%%%%%%%%%%%%%%%
% MEDIAN FLAG 3 PHASE 1
par.mu_max = 0.25266;
par.K_s = 2.6505 ;
par.k_d = 0.02169;
par.Y_xs = 0.4984;
par.Y_xco2 = 0.6047;
end

function par = param3()
%%%%%%%%%%%%%%%
par.S_in = 100;
par.F_in = 0;
par.F_out = 0;
par.q_air = 2;
%%%%%%%%%%%%%%%
% % MEDIAN FLAG 3 PHASE 3
par.mu_max = 2.854E-22;
par.K_s = 20 ;
par.k_d = 0.01258;
par.Y_xs = 9.1688E-22;
par.Y_xco2 = 9.4948E-22;
end

% function par = param2()
% %%%%%%%%%%%%%%%
% par.S_in = 100;
% par.F_in = 0.226825385;
% par.F_out = 0;
% par.q_air = 2;
% %%%%%%%%%%%%%%%
% par.mu_max = 0.310638327;
% par.K_s = 1.56059E-09 ;
% par.k_d = 4.01657E-11;
% par.Y_xs = 1;
% par.Y_xco2 = 0.76799623;
% end

% Function to define all initial conditions
function x0 = init_cond()
% Define initial conditions
Volume_0 = 1.5; %L
Biomass_0 = 1.211; % g/L
Substrate_0 = 20; % g/L
CO2_0 = 0.0282; % percent

% Save all initial conditons in one vector
x0 = [Volume_0; Biomass_0; Substrate_0; CO2_0];
end

function [tX, tS, CO2m, tCO2, S, X] = expdata()
    % function [t, x] = ode45model()
    load measurements
    %using only the timepoints of measurement for substrate and biomass for CO2 since the CO2 measurements are continous
    measured = load("XSCt.mat");
    CO2m = measured.sol(:,3);
    tCO2 = measured.tspan;
end

function [Xexp Sexp Cexp] = expdataphase3()
    load XSCt3 sol
    Xexp = sol(:,1);
    Sexp = sol(:,2);
    Cexp = sol(:,3);
end