clear; clc; close all;

visualization =true;

T = 10;
theta1_initial = 30;
theta2_initial = 45;


y0 = [deg2rad(theta1_initial), deg2rad(theta2_initial) 0, 0];

[t,y] = ode45(@ode_rrbot,[0,T],y0);

K = [    37.2190   -0.3267   14.4845    0.7429;
    9.5519    3.8883    4.1765    0.7474];

U = -K*y';

% tau1 = K(1,1)*theta1 + K(1,2)*theta2 + K(1,3)*theta1_dot + K(1,4)*theta2_dot;
% tau2 = K(2,1)*theta1 + K(2,2)*theta2 + K(2,3)*theta1_dot + K(2,4)*theta2_dot;


%% plots
figure
hold on
subplot(2,2,1)
plot(t,y(:,1))
xlabel('Time step')
ylabel('rad')
title('theta1')

subplot(2,2,2)
plot(t,y(:,2))
xlabel('Time step')
ylabel('rad')
title('theta2')

subplot(2,2,3)
plot(t,y(:,3))
xlabel('Time step')
ylabel('rad/s')
title('theta1-dot')

subplot(2,2,4)
plot(t,y(:,4))
xlabel('Time step')
ylabel('rad/s')
title('theta2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(t,U(1,:))
xlabel('Time step')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(t,U(2,:))
xlabel('Time step')
ylabel('Nm')
title('tau2')
hold off
pause(10)
close all

%% visualization
figure
if(visualization)

    l1 = 1;
    l2 = 1;
    
    x1_pos= l1*sin(y(:,1));
    x2_pos= l1*sin(y(:,1)) + l2*sin(y(:,1)+y(:,2));
    y1_pos= l1*cos(y(:,1));
    y2_pos= l1*cos(y(:,1)) + l2*cos(y(:,1)+y(:,2));
    for i=1:size(y)
        plot([0 x1_pos(i) x2_pos(i)],[0 y1_pos(i) y2_pos(i)], 'LineWidth',4.0)
        hold on
        plot(x1_pos(i),y1_pos(i),'.','MarkerSize',24.0)
        xlim([-2.25 2.25])
        ylim([-2.25 2.25])
        pause(0.000075)
        hold off
    end
end
