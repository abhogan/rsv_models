% Simple SEIR model for RSV transmission - with births and deaths and with
% two age classes, sinusoidally forcing the two transmission coefficients 
% betaA and betaB.
% Author: Alexandra B Hogan
% Year: 2013
% Published in AB Hogan, GN Mercer, K Glass & HC Moore (2013) "Modelling 
% the seasonality of respiratory syncytial virus in young children". 
% Proceedings of the 20th International Congress on
% Modelling and Simulation, Adelaide, Australia 1-6 December 2013.

clf
clear

tend=100*52;   % end time of calculations in weeks
dt=0.1;
tspan=0:dt:tend;  % force solutions on weekly basis

% disease parameters - all in weeks
gamma=1/1.43; % infectious period of 10 days
delta=1/0.57; % latency period of 4 days
nu=1/28.57; % immunity period of 200 days
mu=346/107816; %weekly birth/death rate
eta=1/52; %ageing rate
beta10=3.2; % transmission coefficient in older age group
beta11=0.6; % seasonality in younger age group
beta20=0.75*beta10; % transmission coefficient in younger age group
beta21=beta11; % seasonality in older age group

% initial values
I10=0.0001;
E10=0;
S10=1-I10;
R10=0;
I20=I10;
E20=0;
S20=S10;
R20=0;

% solve ODES using Matlab's ode45 integrator
param=[gamma delta nu mu eta beta10 beta11 beta20 beta21];
[t,y1]=ode45(@SEIRS_ODE_2ages,tspan,[S10 E10 I10 R10 S20 E20 I20 R20],[],param);

% plot Infectives
figure(1)
plotstart = tend-52*12;
plotend = tend;

%plots the sum of infectives in both age classes
inf1 = y1(:,3)./(y1(:,1)+y1(:,2)+y1(:,3)+y1(:,4));
inf2 = y1(:,7)./(y1(:,1)+y1(:,2)+y1(:,3)+y1(:,4));
ymax = max([max(inf1(tend-5*52/dt:tend)) max(inf2(tend-5*52/dt:tend))]) * 1.5;

plot(t-plotstart,inf1,'b-','LineWidth',1.5)%
hold on
plot(t-plotstart,inf2,'Color',[0 0.5 0], 'LineStyle', '--','LineWidth',1.6)
hold on

xlabel('Weeks','FontSize',16)
ylabel('Proportion Infected','FontSize',16)

% note legend definitions not updating dynamically - need to adjust
% manually

h_legend=legend('I_1:  \beta_0=3.2  \beta_1=0.6', 'I_2:  \beta_0=2.4  \beta_1=0.6');
%set(h_legend,'FontSize',16);
axis([0 300 0 ymax])
set(gca, 'xcolor', [0 0 0],'ycolor', [0 0 0],'color', 'none','FontSize',16);
box off
print -djpeg99 SEIRS2betas % create a jpeg of the plot
