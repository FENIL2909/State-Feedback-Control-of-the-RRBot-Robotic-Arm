clear; close; clc;
% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

K = [    37.2190   -0.3267   14.4845    0.7429;
    9.5519    3.8883    4.1765    0.7474];

i = 1;

tic;
t = 0;
while(t < 10)
t = toc;
% read the joint states
jointData = receive(JointStates);
% inspect the "jointData" variable in MATLAB to get familiar with its structure
% design your state feedback controller in the following
Theta1 = wrapTo2Pi(jointData.Position(1));
Theta2 = wrapTo2Pi(jointData.Position(2));
X = [Theta1;Theta2;jointData.Velocity(1);jointData.Velocity(2)];
U = -(K*X);
tau1.Data = U(1);
tau2.Data = U(2);
send(j1_effort,tau1);
send(j2_effort,tau2);
% you can sample data here to be plotted at the end
X1(i) = Theta1;
X2(i) = Theta2;
X3(i) = jointData.Velocity(1);
X4(i) = jointData.Velocity(2);
TAU1(i) = U(1);
TAU2(i) = U(2);

time(i) = t;
i = i+1;

end
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;

%% plots
figure
hold on
subplot(2,2,1)
plot(time,X1)
xlabel('Time step')
ylabel('rad')
title('theta1')

subplot(2,2,2)
plot(time,X2)
xlabel('Time step')
ylabel('rad')
title('theta2')

subplot(2,2,3)
plot(time,X3)
xlabel('Time step')
ylabel('rad/s')
title('theta1-dot')

subplot(2,2,4)
plot(time,X4)
xlabel('Time step')
ylabel('rad/s')
title('theta2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(time,TAU1)
xlabel('Time step')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(time,TAU2)
xlabel('Time step')
ylabel('Nm')
title('tau2')
hold off