% clear
function [t, x] = ode45model()

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
tspan = [21.366701 24.3667];
par.F_in = 0.5/(tspan(2)-tspan(1));
[t2,x2] = ode45(@diff_eq,tspan,x0,options,par);

%% End of additions
x0 = x2(end,:);
par.F_in = 0;
par.F_out = 0;
tspan = [24.366701 50];
[t3,x3] = ode45(@diff_eq,tspan,x0,options,par);

x = [x1; x2; x3];
t = [t1; t2; t3];

Volume = x(:, 1);
Biomass = x(:, 2);
Substrate = x(:, 3);
CO2 = x(:, 4);
bp = x(:,5);
bp2 = x(:,6);

plotting(t,x)

% Function to define your differential equation
function dxdt = diff_eq(t,x,par)

Volume = x(1);
Biomass = x(2);
Substrate = x(3);
CO2 = x(4);
bp = x(5);
bp2 = x(6);

%defining lower boundaries for all states
Volume(Volume<0.00001) = 0.00001;
Biomass(Biomass<0) = 0;
Substrate(Substrate<0) = 0;
CO2(CO2<0) = 0;
bp(bp<0) = 0;
bp2(bp2<0) = 0;

%% Define differential equation
Substrate_dot = (par.F_in./Volume).*(par.S_in-Substrate) - par.mu_maxs.*Substrate./(par.K_s + par.gamma.*bp^par.c + Substrate).*Biomass - par.zeta.*Substrate.*Biomass;

bp2_dot= - bp2* par.beta .* 1/(par.K_s.^par.p + Substrate.^par.p); 
bp_dot = -par.alpha.*bp - bp2_dot; 

% Biomass_dot = -(par.F_in./Volume).*Biomass + S_cons.*Biomass - par.k_d.*Biomass;
Biomass_dot = -(par.F_in./Volume).*Biomass + par.mu_maxx.*(1-bp./(bp+par.K_m1).*par.epsilon).*(Substrate./(Substrate+par.K_m2)).*Biomass -par.k_d.*Biomass;

CO2_dot = par.mu_maxc.*(Substrate./(Substrate+par.K_m2)).*Biomass -par.q_air.*CO2;
Volume_dot = par.F_in;

dxdt = [Volume_dot; Biomass_dot; Substrate_dot; CO2_dot; bp_dot; bp2_dot];
end

% Function to define all your constant parameters
function par = param()		
% Do not change%
par.S_in = 100;%
par.F_in = 0;  %
par.F_out = 0; %
par.q_air = 2; %
%%%%%%%%%%%%%%%%
% 
par.epsilon = 0.4;
par.k_d = 0.02;
par.K_s = 0.1;
par.p = 8;
par.alpha = 0.12;
par.beta = 0.05;
% par.gamma = 10;
% par.c = 0.63;
% par.mu_maxs = 0.226;
% par.mu_maxx = 0.8;
% par.mu_maxc = 1.2143;
% par.K_m1 = 8;
% par.K_m2 = 38;
par.zeta = 0.01;

a = [11.51992826	0.035285763	0.358040072	0.454093958	0.654094624	9.967065343	13.49702579]; %flag 3 MonteCarlo_500

par.gamma = a(1);
par.c = a(2);
par.mu_maxs = a(3);
par.mu_maxx = a(4);
par.mu_maxc = a(5);
par.K_m1 = a(6);
par.K_m2 = a(7);

% par.factor = a(8);
% par.k_d = a(9);
% par.K_s = a(10);
% par.p = a(11);
% par.alpha = a(12);
% par.beta = a(13);
% par.zeta = a(14);

end

% Function to define all initial conditions
function x0 = init_cond()
% Define initial conditions
Volume_0 = 1.5; %L
Biomass_0 = 1.2; % g/L
Substrate_0 = 20; % g/L
CO2_0 = 0; % percent
bp_0 = 0;
bp2_0 = 10;

% Save all initial conditons in one vector
x0 = [Volume_0; Biomass_0; Substrate_0; CO2_0; bp_0; bp2_0];
end

end
function plotting(t,x)
    load measurements
    %using only the timepoints of measurement for substrate and biomass for CO2 since the CO2 measurements are continous
    measured = load("XSCt.mat");
    CO2m = measured.sol(:,3);
    tCO2 = measured.tspan;

    % Volume = x(:, 1);
    Biomass = x(:, 2);
    Substrate = x(:, 3);
    CO2 = x(:, 4);
    bp = x(:,5);
    bp2= x(:,6);

    % Plot figure
    figure(1)
    plot(t, Biomass,"Color","blue")
    hold on
    plot(t, Substrate,"Color","red")
    plot(t, CO2, "Color", "#77AC30")
    plot(tX,X,"o","Color","blue")
    plot(tS,S,"o","Color","red")
    plot(tCO2,CO2m,"o","Color", "#77AC30")
    hold off
    ylim([0 25])
    xlabel("Time [h]")
    legend("X [g/L]", "S [g/L]", "CO_2 [%]", "Measured X [g/L]", "Measured S [g/L]", "Measured CO_2 [%]")
    % saveas(gcf,path,"png")
    % close(1)
    
    figure(2)
    plot(t, bp, "Color", "black")
    hold on
    plot(t,bp2,'--k')
    hold off
    ylim([0 10])
    xlabel("Time [h]")
    legend("bp [g/L] ","bp_2 [g/L] ")
end