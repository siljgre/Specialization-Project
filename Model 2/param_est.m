function [results, EXITFLAG] = param_est(initParams, loBound, upBound) 
% Estimation of the parameters by non-linear least squares optimization
% problem. 

sol=[0 0 0 0]; 
tspan = 0;

load XSCt sol tspan
exp=sol;

%får ikkje den feilen med default metode, men resultat ser ikkje bra ut
% opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
opts = optimset('TolFun', 1e-12, 'TolX', 1e-6, 'MaxIter', 150,'Algorithm', 'levenberg-marquardt', 'Diagnostics', 'off', 'Display', 'iter'); 

%minimize residual with lsqnonlin

[res, RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);

function R = residual(b)
    a = 10.^b;

    [T ,Y]=reactionsolve (a);

    % Y = V X S C = solution of ode15s
    % exp = X S C bp = experimental values
    residual = exp - Y(:,2:4);
    R = residual(:);
    results=a;

    % Plot results
    plotting(T,Y,tspan,exp)
end

function [ T, Y ] = reactionsolve( a )

% parameters to estimate
par.c = a(1); % if less than 1 then the term par.gamma.*bp^par.c becomes large. 
par.p = a(2);
par.alpha = a(3);
par.beta = a(4);
par.gamma = a(5);
par.zeta = a(6);
par.K_m = a(7);

% Initial conditions
x0 = init_cond();

%% Call ode15s
% Phase 1
addition = 0;
tspan0=tspan(1:17);
[t1,y1] = ode15s(@reaction, tspan0, x0, [], addition);
x0 = y1(end,:);
t1(17,:) = [];
y1(17,:) = [];

% Phase 2
addition = 1; 
tspan1 = tspan(17:19);
[t2,y2] = ode15s(@reaction, tspan1, x0, [], addition); %at error: [t2, y1] is only one element (t2=21.5)
x0 = y2(end,:); 
t2(end,:) = []; %at error: then t2 = []
y2(end,:) = [];

% Phase 3
addition = 0;
tspan2 = tspan(19:end);
[t3,y3] = ode15s(@reaction, tspan2, x0, [], addition);
T=[t1;t2;t3];
Y=[y1;y2;y3];

function dx = reaction(t, x, addition)
%%%%%%%%%%%%%%%%
% Do not change%
par.S_in = 100;
par.F_in = 0;  
par.F_out = 0; 
par.q_air = 2; 
%%%%%%%%%%%%%%%%
if addition
    par.F_in=0.5/2.8667; %0.5 L divided by the timespan of addition (tspan(19)-tspan(17))
end

par.mu_max = 0.27 *0.6;
par.K_s = 2.5;
par.k_d = 0.018;
par.Y_xs = 0.42042; 
par.Y_xco2 = 0.604716013;
par.c2 = 6;
par.c3 = 0.4/0.6;

% par.c = 0.1;
% par.p = 6.5;
% par.alpha = 0.1;
% par.beta = 2.3386;
% par.gamma = 57;
% par.zeta = 0.01; 
% par.K_m = 0.01;

%% Unpacking states from array
Volume = x(1);
Biomass = x(2);
Substrate = x(3);
CO2 = x(4);
bp = x(5);
bp2 = x(6);
alpha = x(7);
alpha2 = x(8);

%% No states should go below zero
Volume(Volume<0.00001) = 0.00001; %skjer hvis inni parantes er riktig
Biomass(Biomass<0) = 0;
Substrate(Substrate<0) = 0;
CO2(CO2<0) = 0;
bp(bp<0)=0;
bp2(bp2<0)=0;
alpha(alpha<0)=0;
alpha2(alpha<0)=0;

% % Stop emptying the tank when the tank is empty, to ensure that volume does
% % not reach negative values
% if Volume == 0.00001
%     par.F_out = 0;
% end

%% Model equations
S_cons = par.mu_max*Substrate./(par.K_s + 1*alpha*par.gamma.*bp^par.c + Substrate); %factor for S consumption for biomass growth
S_consX= alpha*S_cons*par.c2 + alpha2*par.c3*par.mu_max.*Substrate./(par.K_s + Substrate)./par.Y_xs; %factor for S consumption for biomass growth
Substrate_dot = (par.F_in./Volume).*(par.S_in-Substrate) - S_cons.*(Biomass./(par.Y_xs)) - par.zeta.*Substrate.*Biomass;

Biomass_dot = -(par.F_in./Volume).*Biomass + S_consX.*Biomass - par.k_d.*Biomass;

CO2_dot = S_consX.*(Biomass./par.Y_xco2)-par.q_air.*CO2;
Volume_dot = par.F_in;

bp_dot = -par.alpha.*bp + bp2* par.beta .* 1/(par.K_m.^par.p + Substrate.^par.p);
bp2_dot= - bp2* par.beta .* 1/(par.K_m.^par.p + Substrate.^par.p); % an "intermediate" between bp and bp2?

alpha_dot= alpha2.* 1/(par.K_m.^par.p + Substrate.^par.p); %alpha is 1 during starvation, 0 otherwise (gradually 1 and zero)
alpha2_dot = -alpha2.* 1/(par.K_m.^par.p + Substrate.^par.p); %alpha2 is 0 during starvation, 1 otherwise

dx = zeros(8,1);
dx(1)=Volume_dot;
dx(2)=Biomass_dot;
dx(3)=Substrate_dot;
dx(4)=CO2_dot;
dx(5)=bp_dot;
dx(6)=bp2_dot;
dx(7)=alpha_dot;
dx(8)=alpha2_dot;
end
end
end

function x0 = init_cond()
% Function to define all initial conditions
% Define initial conditions
Volume_0 = 1.5; %L
Biomass_0 = 1.2; % g/L
Substrate_0 = 20; % g/L
CO2_0 = 0; % percent
bp_0 = 0;
bp2_0=5;
alpha_0=0;
alpha2_0=1;

% Save all initial conditons in one vector
x0 = [Volume_0; Biomass_0; Substrate_0; CO2_0; bp_0; bp2_0;alpha_0;alpha2_0];
end

function plotting(T,Y,tspan,exp) 
    phase1 = 1:17; %index. 1 is first measurement, 17 is last measurement before feeding
    phase2 = 17:19;
    phase3 = 19:29;
    
    %% Biomass
    %phase1
    subplot(4,1,1);
    plot(T(phase1),Y(phase1,2));
    hold on
    plot(tspan(phase1),exp(phase1,1),'o');
    %phase2
    plot(T(phase2),Y(phase2,2));
    plot(tspan(phase2),exp(phase2,1),'o');
    %phase3
    plot(T(phase3),Y(phase3,2));
    plot(tspan(phase3),exp(phase3,1),'o');
    xlabel('Time');
    ylabel('Biomass');
    hold off
    % str = sprin2wstf('%g | ', a);
    % title(str);

    %% Substrate
    %phase1
    subplot(4,1,2);
    plot(T(phase1),Y(phase1,3));
    hold on
    plot(tspan(phase1,1),exp(phase1,2),'o');
    %phase2
    plot(T(phase2),Y(phase2,3));
    plot(tspan(phase2,1),exp(phase2,2),'o');
    %phase3
    plot(T(phase3),Y(phase3,3));
    plot(tspan(phase3,1),exp(phase3,2),'o');
    hold off
    xlabel('Time');
    ylabel('Substrate');
    hold off;

    %% CO2
    subplot(4,1,3);
    plot(T(1:17),Y(1:17,4));
    hold on;
    plot(tspan(1:17,1),exp(1:17,3),'o');
    %phase2
    plot(T(17:19),Y(17:19,4));
    plot(tspan(17:19,1),exp(17:19,3),'o');
    %phase3
    plot(T(phase3),Y(phase3,4));
    plot(tspan(phase3,1),exp(phase3,3),'o');
    hold off
    xlabel('Time');
    ylabel('CO2');
    hold off;

    %% bp
    subplot(4,1,4);
    plot(T(phase1),Y(phase1,5));
    hold on
    %phase2
    plot(T(phase2),Y(phase2,5));
    %phase3
    plot(T(phase3),Y(phase3,5));
    xlabel('Time');
    ylabel('By-products');
    hold off;
    drawnow;
end
