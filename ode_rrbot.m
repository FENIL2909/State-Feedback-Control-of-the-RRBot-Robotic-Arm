%% ode_rrbot

function dX = ode_rrbot(t,X)

M1 = 1; %Kg
M2 = 1; %Kg
L1 = 1; %m
L2 = 1; %m
r1 = 0.45; %m
r2 = 0.45; %m
I1 = 0.084; %Kg.m2
I2 = 0.084; %Kg.m2
g = 9.81; %m/s2

dX = zeros(4,1);
X = num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});
% theta1 = wrapTo2Pi(theta1);
% theta2 = wrapTo2Pi(theta2);

K = [   37.2190   -0.3267   14.4845    0.7429;
    9.5519    3.8883    4.1765    0.7474];

U = -K*[theta1; theta2; theta1_dot; theta2_dot];
tau1 = U(1);
tau2 = U(2);

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*tau1 - I2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 + L1*M2^2*g*r2^2*sin(theta1) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) - L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - L1^2*M2*tau2 - M1*r1^2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 - L1^2*M2^2*g*r2*sin(theta1 + theta2) + L1*M2^2*g*r2^2*sin(theta1) - I1*M2*g*r2*sin(theta1 + theta2) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) + L1*M2*r2*tau1*cos(theta2) - 2*L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1^3*M2^2*r2*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + L1^2*M2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + L1^2*M2^2*g*r2*cos(theta2)*sin(theta1) - M1*M2*g*r1^2*r2*sin(theta1 + theta2) + I1*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + L1*M1*M2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2) + L1*M1*M2*g*r1*r2*cos(theta2)*sin(theta1))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
end