function dy = SEIRS_ODE(t,y,param)

% ode function for the SEIR model - with births and deaths
% Author: Alexandra B Hogan (2013)

S = y(1); %susceptible
E = y(2); %exposed
I = y(3); %infected
R = y(4); %recovered

gamma = param(1); %1/gamma is av recovery period
delta = param(2); %1/delta is average latency period
nu = param(3); %1/nu is average duration of immunity
mu = param(4); %birth rate/ageing rate
beta0 = param(5);
beta1 = param(6);

% beta varies sinusoidally
beta=beta0*(1+beta1*sin(2*pi*t/52));

dy=[mu - beta*S*I - mu*S + nu*R;
    beta*S*I - delta*E - mu*E;
    delta*E - mu*I - gamma*I;
    gamma*I - mu*R - nu*R];
            