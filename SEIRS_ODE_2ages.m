function dy = SEIRS_ODE_2ages(t,y,param)

% ode function for the SEIR model for RSV - with births and deaths and 
% with two age classes and two transmission parameters
% Author: Alexandra B Hogan (2013)

S1 = y(1); %susceptible
E1 = y(2); %exposed
I1 = y(3); %infected
R1 = y(4); %recovered
S2 = y(5);
E2 = y(6);
I2 = y(7);
R2 = y(8);

gamma = param(1); %1/gamma is av recovery period
delta = param(2); %1/delta is average latency period
nu = param(3); %1/nu is average duration of immunity
mu = param(4); %birth rate
eta = param(5); %ageing rate
beta10 = param(6);
beta11 = param(7);
beta20 = param(8);
beta21 = param(9);

% beta varies sinusoidally
beta1=beta10*(1+beta11*sin(2*pi*t/52));
beta2=beta20*(1+beta21*sin(2*pi*t/52));
N1=S1+E1+I1+R1;
N2=S2+E2+I2+R2;
dy=[mu - beta1*S1*(I1 + I2) - eta*S1 + nu*R1
    beta1*S1*(I1 + I2) - E1*(delta + eta)
    delta*E1 - I1*(eta + gamma)
    gamma*I1 - R1*(eta + nu)
    eta*S1 - S2*beta2*(I1 + I2) - eta*S2 + nu*R2
    eta*E1 + S2*beta2*(I1 + I2) - E2*(delta + eta)
    eta*I1 + delta*E2 - I2*(eta + gamma)
    eta*R1 + gamma*I2 - R2*(eta + nu)];
