% ode45 method

function [t, x] = ode45model()

load measurements
%using only the timepoints of measurement for substrate and biomass for CO2 since the CO2 measurements are continous
measured = load("XSCt.mat");
CO2m = measured.sol(:,3);
tCO2 = measured.tspan;

% Initialize initial conditions
x0 = init_cond();
% Initialize parameters
par = param();
% Define timespan differential equations are simulated for
tspan = [0 21.3667];
% Initialize options (ALWAYS EMPTY!)
options = [];

% Run numerical solver
[t1,x1] = ode45(@diff_eq,tspan,x0,options,par);


%% Addition of substrate
x0 = x1(end,:);
% par.F_in = 0.5/(25-21.3667);
% tspan = [21.366701 24.3667];
tspan = [21.366701 24.3667];
par.F_in = 0.5/(tspan(2)-tspan(1));
[t2,x2] = ode45(@diff_eq,tspan,x0,options,par);

%% After addition
x0 = x2(end,:);
par.F_in = 0;
tspan = [24.366701 50];
[t3,x3] = ode45(@diff_eq,tspan,x0,options,par);

% Concatenation of arrays
x = [x1; x2; x3];
t = [t1; t2; t3];

Volume = x(:, 1);
Biomass = x(:, 2);
Substrate = x(:, 3);
CO2 = x(:, 4);

% Plot figure
% figure()
% plot(t, Volume)
% hold on
% plot(t, Biomass)
% plot(t, Substrate)
% plot(t, CO2)
% hold off
plot(t, Biomass,"Color","blue")
hold on
plot(t, Substrate,"Color","red")
plot(t, CO2, "Color", "#77AC30")
plot(tX,X,"o","Color","blue")
plot(tS,S,"o","Color","red")
plot(tCO2,CO2m,"o","Color", "#77AC30")
hold off
xlabel("Time (h)")
ylim([0 25])
legend("Biomass", "Glucose", "CO2")
legend("Biomass [g/L]", "Glucose [g/L]", "CO2 [%]")
end

% Function to define your differential equation
function dxdt = diff_eq(t,x,par)
Volume = x(1);
Biomass = x(2);
Substrate = x(3);
CO2 = x(4);

Volume(Volume<0.00001) = 0.00001; %skjer hvis inni parantes er riktig
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
function par = param()
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

% Function to define all initial conditions
function x0 = init_cond()
% Define initial conditions
Volume_0 = 1.5; %L
Biomass_0 = 1.2; % g/L
Substrate_0 = 20; % g/L
CO2_0 = 0; % percent

% Save all initial conditons in one vector
x0 = [Volume_0; Biomass_0; Substrate_0; CO2_0];
end