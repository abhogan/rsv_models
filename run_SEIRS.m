% Simple SEIRS model for RSV transmission - with births and deaths, for one
% age class
% Author: Alexandra B Hogan
% Year: 2013
% Published in AB Hogan, GN Mercer, K Glass & HC Moore (2013) "Modelling 
% the seasonality of respiratory syncytial virus in young children". 
% Proceedings of the 20th International Congress on
% Modelling and Simulation, Adelaide, Australia 1-6 December 2013.

clear all
clf
tend=52*80;   % end time of calculations in weeks
tspan=0:0.1:tend;  % force solutions on weekly basis

% disease parameters
gamma=1/1.4;
delta=1/0.57;
nu=1/28.57;
mu=346/107816; %weekly birth rate
beta0=1.1;
beta1=0.6;

% initial values
I0=0.001;
S0=1-I0;
E0=0;
R0=0;

% solve ODES using Matlab's ode45 integrator
param=[gamma delta nu mu beta0 beta1];
[t,y1]=ode45(@SEIRS_ODE,tspan,[S0 E0 I0 R0],[],param);

% plot Infectives
figure(1)
plotstart=1100;
plotend=tend;
plot(t-plotstart,y1(:,3),'b-','LineWidth',1.5)
hold on
xlabel('Weeks','FontSize',16)
ylabel('Proportion Infected','FontSize',16)
axis([0 350 0 0.2])

figure(2)
plot(t,y1(:,3))
print -djpeg99 SIRexampleI