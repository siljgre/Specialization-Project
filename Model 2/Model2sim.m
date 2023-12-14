
function [t,x] = Model2sim()
% Initialize initial conditions and parameters
x0 = init_cond();
par = param();
% Initialize options (ALWAYS EMPTY!)
options = [];

%% Phase 1 (Batch)
tspan = [0 21.3667];
[t1,x1] = ode15s(@diff_eq,tspan,x0,options,par);

%% Phase 2 (Feeding)
x0 = x1(end,:);
tspan = [21.366701 24.3667];
%adding 0.5L of feed over time interval tspan
par.F_in = 0.5/(tspan(2)-tspan(1));
[t2,x2] = ode15s(@diff_eq,tspan,x0,options,par);

%% Phase 3 (After feeding)
x0 = x2(end,:);
par.F_in = 0;
tspan = [24.366701 50];
[t3,x3] = ode15s(@diff_eq,tspan,x0,options,par);

% Saving solutions to one state array x and one time array t
x = [x1; x2; x3];
t = [t1; t2; t3];

%% Plot solutions
plotting(t,x)
end

%% Functions
function dxdt = diff_eq(t,x,par)
Volume = x(1);
Biomass = x(2);
Substrate = x(3);
CO2 = x(4);
bp = x(5);
bp2=x(6);
alpha=x(7);
alpha2=x(8);

%defining lower boundaries for all states
Volume(Volume<0.00001) = 0.00001; %skjer hvis inni parantes er riktig
Biomass(Biomass<0) = 0;
Substrate(Substrate<0) = 0;
CO2(CO2<0) = 0;
bp(bp<0) = 0;

%bp2(bp2<0.00002)=0.00002;
%% Define differential equation
S_cons = par.mu_max*Substrate./(par.K_s + 1*alpha*par.gamma.*bp^par.c + Substrate); %factor for S consumption for biomass growth
S_consX= (alpha*S_cons*par.c2) + alpha2*par.c3*par.mu_max.*Substrate./(par.K_s + Substrate)./par.Y_xs; %factor for S consumption for biomass growth
Substrate_dot = (par.F_in./Volume).*(par.S_in-Substrate) - S_cons.*(Biomass./(par.Y_xs)) - par.zeta.*Substrate.*Biomass;

Biomass_dot = -(par.F_in./Volume).*Biomass + S_consX.*Biomass - par.k_d.*Biomass;

CO2_dot = S_consX.*(Biomass./par.Y_xco2)-par.q_air.*CO2;

Volume_dot = par.F_in;

bp_dot = -par.alpha.*bp + bp2* par.beta .* 1/(par.K_m.^par.p + Substrate.^par.p);
bp2_dot= - bp2* par.beta .* 1/(par.K_m.^par.p + Substrate.^par.p); % an "intermediate" between bp and bp2?

alpha_dot= alpha2.* 1/(par.K_m.^par.p + Substrate.^par.p); %alpha is 1 during starvation, 0 otherwise (gradually 1 and zero)
alpha2_dot = -alpha2.* 1/(par.K_m.^par.p + Substrate.^par.p); %alpha2 is 0 during starvation, 1 otherwise

% Save differential equation in output vector
dxdt = [Volume_dot; Biomass_dot; Substrate_dot; CO2_dot; bp_dot;bp2_dot;alpha_dot;alpha2_dot];
end

% Function to define all your constant parameters
function par = param()		
% % Do not change%
par.S_in = 100;%
par.F_in = 0;  %
par.F_out = 0; %
par.q_air = 2; % 
%%%%%%%%%%%%%%%%
par.mu_max = 0.27 *0.605;
par.K_s = 2.5;
par.k_d = 0.018;
par.Y_xs = 0.42042;
par.Y_xco2 = 0.604716013;
par.c2 = 5.8;
par.c3 = 0.4/0.6;
% par.c = 0.03;
% par.p = 4;
% par.alpha = 0.1;
% par.beta = 2.3386;
% par.gamma = 57;
% par.zeta = 0.01; 
% par.K_m = 0.01;

% a = [3.59906E-05	2.221735648	0.073461455	4.492547094	44.24452813 0.018260317	1]; %median AnalysisMod2 (100 iterations)
a = [0.000173692	2.221446836	0.082733776	4.157078656	44.24201569	0.018260438	0.999999966]; %median MonteCarlo_500

par.c = a(1);
par.p = a(2);
par.alpha = a(3);
par.beta = a(4);
par.gamma = a(5);
par.zeta = a(6);
par.K_m = a(7);
end

% Function to define all initial conditions
function x0 = init_cond()
% Define initial conditions
Volume_0 = 1.5; %L
Biomass_0 = 1.2; % g/L
Substrate_0 = 20; % g/L
CO2_0 = 0; % percent
bp_0 = 0;
bp2_0=10;
alpha_0=0;
alpha2_0=1;

% Save all initial conditons in one vector
x0 = [Volume_0; Biomass_0; Substrate_0; CO2_0; bp_0; bp2_0;alpha_0;alpha2_0];
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
    alpha=x(:,7);
    alpha2=x(:,8);

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
    xlabel("Time (h)")
    legend("Biomass [g/L]", "Substrate [g/L]", "CO_2 [%]")
    % saveas(gcf,path,"png")
    % close(1)
    
    figure(2)
    subplot(2,1,1)
    plot(t, Biomass,"Color","blue")
    hold on
    plot(t, Substrate,"Color","red")
    plot(t, CO2, "Color", "#77AC30")
    plot(tX,X,"o","Color","blue")
    plot(tS,S,"o","Color","red")
    plot(tCO2,CO2m,"o","Color", "#77AC30")
    hold off
    xlabel("Time (h)")
    legend("Biomass [g/L]", "Substrate [g/L]", "CO_2 [%]")
    
    subplot(2,1,2)
    plot(t, bp, "Color", "black")
    hold on
    plot(t,bp2,'--k')
    plot(t,alpha,'-b')
    plot(t,alpha2,"--b")
    hold off
    legend("bp","bp_2","{\lambda_1}", "{\lambda_2}") %used lambda in report so alpha = lambda
end