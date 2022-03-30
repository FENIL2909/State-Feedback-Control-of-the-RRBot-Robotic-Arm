%% Dynamic Equations
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot tau1 tau2 g 'real'
syms M1 M2 I1 I2 L1 L2 r1 r2 'real'
syms x1_dot y1_dot x2_dot y2_dot 'real'



x1_dot = (theta1_dot)*(r1)*(cos(theta1));
y1_dot = -(theta1_dot)*(r1)*(sin(theta1));

x2_dot = (theta1_dot)*(L1)*(cos(theta1)) + (theta1_dot + theta2_dot)*(r2)*(cos(theta1 + theta2));
y2_dot = -(theta1_dot)*(L1)*(sin(theta1)) - (theta1_dot + theta2_dot)*(r2)*(sin(theta1 + theta2));

K1 = (1/2)*(I1)*(theta1_dot*theta1_dot) + (1/2)*(M1)*((x1_dot*x1_dot) + (y1_dot*y1_dot));
K2 = (1/2)*(I2)*((theta2_dot + theta1_dot)*(theta2_dot + theta1_dot)) + (1/2)*(M2)*((x2_dot*x2_dot) + (y2_dot*y2_dot));

P1 = M1*g*r1*cos(theta1);
P2 = M2*g*(L1*(cos(theta1)) + r2*(cos(theta1 + theta2)));
L = K1 + K2 - P1 - P2;

u = [tau1;tau2];
q = [theta1;theta2];
dq = [theta1_dot; theta2_dot];
ddq = [theta1_ddot; theta2_ddot];

DL_Dq = gradient(L,q);  % used gradient instead of jacobian to keep matrix size consistent
DL_Ddq = gradient(L,dq); % used gradient instead of jacobian to keep matrix size consistent
dDL_dtDdq = jacobian(DL_Ddq,[q;dq])*[dq;ddq];

EOM = simplify(dDL_dtDdq - DL_Dq -u);

%% Finding Equilibrium Points

EOM_t = subs(EOM,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, tau1, tau2],[0,0,0,0,0,0]);
%display(EOM_t);
solution1 = solve(EOM_t ==0, [theta1,theta2]);

fprintf("********** EQUILIBRIUM POINTS **********\n");

fprintf("values of theta1:\n")
disp(solution1.theta1);

fprintf("values of theata2:\n")
disp(solution1.theta2);

% [0,0], [pi,0], [0,pi]

%% System parameters

M1 = 1; %Kg
M2 = 1; %Kg
L1 = 1; %m
L2 = 1; %m
r1 = 0.45; %m
r2 = 0.45; %m
I1 = 0.084; %Kg.m2
I2 = 0.084; %Kg.m2
g = 9.81; %m/s2


%% State Space Representation
X = sym('X', [4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot;
X(4) = theta2_dot;

eq1 = EOM(1);
eq2 = EOM(2);
solution2 = solve([eq1==0, eq2==0],[theta1_ddot,theta2_ddot]);

% display(solution2.theta1_ddot);
% display(solution2.theta2_ddot);

X1_dot = theta1_dot;
X2_dot = theta2_dot;
X3_dot = (I2*tau1 - I2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 + L1*M2^2*g*r2^2*sin(theta1) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) - L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
X4_dot = -(I2*tau1 - I1*tau2 - I2*tau2 - L1^2*M2*tau2 - M1*r1^2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 - L1^2*M2^2*g*r2*sin(theta1 + theta2) + L1*M2^2*g*r2^2*sin(theta1) - I1*M2*g*r2*sin(theta1 + theta2) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) + L1*M2*r2*tau1*cos(theta2) - 2*L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1^3*M2^2*r2*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + L1^2*M2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + L1^2*M2^2*g*r2*cos(theta2)*sin(theta1) - M1*M2*g*r1^2*r2*sin(theta1 + theta2) + I1*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + L1*M1*M2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2) + L1*M1*M2*g*r1*r2*cos(theta2)*sin(theta1))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);

X_dot = [X1_dot;X2_dot;X3_dot;X4_dot];

fprintf("********** STATE SPACE REPRESENTATION **********\n")
fprintf("X1_dot:")
fprintf("theta1_dot\n")

fprintf("X2_dot:")
fprintf("theta2_dot\n")

fprintf("X3_dot:")
disp(solution2.theta1_ddot)

fprintf("X4_dot:")
disp(solution2.theta2_ddot)

%Linearization
A = jacobian(X_dot,X);
B = jacobian(X_dot,u);

fprintf("********** LINEARIZATION AROUND EQUILIBRIUM POINTS **********\n")



%% find A and B at all equilibrium points
% 1 - [0,0,0,0]
A1 = subs(A,[X(1),X(2),X(3),X(4)],[0,0,0,0]);
B1 = subs(B,[X(1),X(2),X(3),X(4)],[0,0,0,0]);
A1 = double(A1);
B1 = double(B1);
fprintf("----Equilibrium Point (0,0,0,0):\n")
fprintf("State Matrix (A): \n")
disp(A1)
fprintf("Input Matrix (B): \n")
disp(B1)

% 2 - [pi,0,0,0]
A2 = subs(A,[X(1),X(2),X(3),X(4)],[pi,0,0,0]);
B2 = subs(B,[X(1),X(2),X(3),X(4)],[pi,0,0,0]);
A2 = double(A2);
B2 = double(B2);
fprintf("----Equilibrium Point (pi,0,0,0):\n")
fprintf("State Matrix (A): \n")
disp(A2)
fprintf("Input Matrix (B): \n")
disp(B2)


% 3 - [0,pi,0,0]
A3 = subs(A,[X(1),X(2),X(3),X(4)],[0,pi,0,0]);
B3 = subs(B,[X(1),X(2),X(3),X(4)],[0,pi,0,0]);
A3 = double(A3);
B3 = double(B3);
fprintf("----Equilibrium Point (0,pi,0,0):\n")
fprintf("State Matrix (A): \n")
disp(A3)
fprintf("Input Matrix (B): \n")
disp(B3)

% 4 - [pi,pi,0,0]
A4 = subs(A,[X(1),X(2),X(3),X(4)],[pi,pi,0,0]);
B4 = subs(B,[X(1),X(2),X(3),X(4)],[pi,pi,0,0]);
A4 = double(A4);
B4 = double(B4);
fprintf("----Equilibrium Point (pi,pi,0,0) [EXTRA SOLUTION - NOT RETURNED BY MATLAB]:\n")
fprintf("State Matrix (A): \n")
disp(A4)
fprintf("Input Matrix (B): \n")
disp(B4)

%% Checking stability at each equilibrium points

fprintf("********** CHECKING STABILITY **********\n")
fprintf("-- Eigen values of (A) around equilibrium point (0,0,0,0):\n")
eigA1 = eig(A1);
disp(eigA1)

fprintf("-- Eigen values of (A) around equilibrium point (pi,0,0,0):\n")
eigA2 = eig(A2);
disp(eigA2)

fprintf("-- Eigen values of (A) around equilibrium point (0,pi,0,0):\n")
eigA3 = eig(A3);
disp(eigA3)

fprintf("-- Eigen values of (A) around equilibrium point (pi,pi,0,0) -- EXTRA SOLUTION:\n")
eigA4 = eig(A4);
disp(eigA4)
%% Checking Controllability at the upward equilibrium point [0,0,0,0]

fprintf("********** INVESTIGATING CONTROLLABILITY AT UPWARD CONFIGURATION **********\n")
fprintf("rank of controllability matrix:")
rankc = rank(ctrb(A1,B1));   % = 4  (controllable)
disp(rankc)

%% State feedback control for upward configuration
syms k1 k2 k3 k4 k5 k6 k7 k8 lambda

% X belong to R^n; here n = 4
% U belong to R^m, here m = 2
% K belong to R^m*n, here, 2*4     [k1,k2,k3,k4;k5,k6,k7,k8]

% control law U = -KX
% Acl = A-Bk

lambda = [-1.4,-2.6,-3.6,-6.7];    %[-2,-4,-3+1i,-3-1i]; [-9,-12,-15,-17]; [-2,-3,-4,-7]; [-1.4,-2.6,-3.6,-6.7];

fprintf("********** K Gains **********")
K = place(A1,B1,lambda)




