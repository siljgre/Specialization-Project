function [results, EXITFLAG] = param_est(initParams,loBound, upBound); 
% Estimation of the parameters by non-linear least squares optimization
% problem. 

% run for example: model_parameter_estimation_tests3([0.1;0.1;0.1]) for 3 parameters. 
% run for example: model_parameter_estimation_tests6([0.1;0.1;0.1; 0.5;0.5;0.05 ]) for 6 parameters. 

% trying to estimate mu_max
sol=[0 0 0]; tspan = 0;
load SXt.mat %only estimating until second addition of sugar
sol = sol(1:16,:);
tspan = tspan(1:16,:);

exp=sol;
%exp is [X S CO2]
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter'); %original
% opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150,'Algorithm', 'levenberg-marquardt', 'Diagnostics', 'off', 'Display', 'iter'); %not original
% opts = optimset('TolFun', 1e-8, 'TolX', 1e-12, 'MaxIter', 100, 'Diagnostics', 'off', 'Display', 'iter');

% The original values are:
% par.mu_max = 0.19445;
% par.K_s = 0.007;
% par.k_d = 0.006;
% par.Y_xs = 0.42042;
% par.Y_xco2 = 0.54308;

% positive feedback reaction from 
% estimate one parameter, mu_max
% loBound = [0]; 
% upBound = [1];

[res, RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);
   function R = residual(b)
        
      a = 10.^b;
         
      [T ,Y ]= reactionsolve (a);
    
      %in this case: looking only at S and X

       residual = exp - Y(:,2:end); %extracting first column (not looking at V)
       R = residual(:);
       results=a;
       
          subplot(3,1,1);
          plot(T,Y(:,3)); %Y(all rows, third column = substrate)
          hold on;
          plot(tspan(:,1),exp(:,2),'o'); %2nd column of exp is substrate
          xlabel('Time');
          ylabel('Substrate');
          hold off;
          str = sprintf('%g | ', a);
          title(str);

          subplot(3,1,2);
          plot(T,Y(:,2)); %Y(all rows, second column = biomass)
          hold on;
          plot(tspan(:,1),exp(:,1),'o'); %1st column of exp is biomass
          xlabel('Time');
          ylabel('Biomass');
          hold off;

          subplot(3,1,3);
          plot(T,Y(:,4));
          hold on;
          plot(tspan(:,1),exp(:,3),'o');
          xlabel('Time');
          ylabel('CO2');
          hold off;

          hold off;
          drawnow;
   end

function [ T, Y ] = reactionsolve( a )


par.mu_max = a(1);
par.K_s = a(2);
par.k_d = a(3);
par.Y_xs = a(4);
par.Y_xco2 = a(5);

x0 = [1.5; 1.2; 20; 0];

[T,Y] = ode45(@reaction, tspan, x0, []);

function dx = reaction(t, x)

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

% Fixed parameters:
% par.mu_max = 0.19445;
% par.K_s = 0.007;
% par.k_d = 0.006;
% par.Y_xs = 0.42042;
% par.Y_xco2 = 0.54308;
par.S_in = 100;
par.F_in = 0;
par.F_out = 0;
par.q_air = 2;

%Diff lign
Volume_dot = par.F_in;
Biomass_dot = -(par.F_in./Volume).*Biomass+par.mu_max.*(Substrate./(par.K_s+Substrate)).*Biomass-par.k_d.*Biomass;
Substrate_dot = (par.F_in./Volume).*(par.S_in-Substrate)-(par.mu_max).*Substrate./((par.K_s+Substrate)).*(Biomass./par.Y_xs);
CO2_dot = par.mu_max.*(Substrate./(par.K_s+Substrate)).*(Biomass./par.Y_xco2)-par.q_air.*CO2;

dx = zeros(4,1);
dx(1)=Volume_dot;
dx(2)=Biomass_dot;
dx(3)=Substrate_dot;
dx(4)=CO2_dot;
    
end


end

end


