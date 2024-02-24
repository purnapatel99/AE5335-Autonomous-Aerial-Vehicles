clear
clc

load asg3_prob4_coefficients.mat

m = 13.5; %kg
g = 9.81; %m/s^2
rho = 1.225; % kg/m^3
S = 0.55; %m^2
CL_0 = 0.28;
CD_0 = 0.03;
CL_alp = 3.45;
CD_alp = 0.30;

t_plot = [];
p_plot = [];
V_pot = [];
chi_plot = [];
gamma_plot = [];
L_plot = [];
D_plot = [];
T_plot = [];
alpha_plot = [];
bank_angle_plot = [];
n = 1;
for t = 0 : 0.1 : 60

   
    p = coeff_p(:,5) + coeff_p(:,4)*t  + coeff_p(:,3)*t^2 + coeff_p(:,2)*t^3 + coeff_p(:,1)*t^4;
    p_dot = coeff_p(:,4)  + 2*coeff_p(:,3)*t + 3*coeff_p(:,2)*t^2 + 4*coeff_p(:,1)*t^3;
    p_2dot = 2*coeff_p(:,3) + 6*coeff_p(:,2)*t + 12*coeff_p(:,1)*t^2;

    V = norm(p_dot);

    chi = atan2(p_dot(2), p_dot(1));

    gamma = asin(-p_dot(3)/V);

    V_dot = (p_dot(1)*p_2dot(1) + p_dot(2)*p_2dot(2) + p_dot(3)*p_2dot(3))/V;

    gamma_dot = (V_dot*p_dot(3) - p_2dot(3)*V)/(V^2*cos(gamma));

    chi_dot = (cos(chi))^2*(p_dot(1)*p_2dot(2) - p_2dot(1)*p_dot(2))/p_dot(1)^2;

    bank_angle = atan2(m*V*chi_dot*cos(gamma), m*V*gamma_dot + m*g*cos(gamma));

    L = m*V*chi_dot*cos(gamma)/sin(bank_angle);

    CL = 2*L/(rho*V^2*S);

    alpha = (CL - CL_0)/CL_alp;

    D = 0.5*rho*V^2*S*(CD_0 + alpha*CD_alp);

    T = m*V_dot + D + m*g*sin(gamma);

    t_plot(:,n) = t;
    p_plot(:,n) = p;
    V_pot(:,n) = V;
    chi_plot(:,n) = chi;
    gamma_plot(:,n) = gamma;
    L_plot(:,n) = L;
    D_plot(:,n) = D;
    T_plot(:,n) = T;
    alpha_plot(:,n) = alpha;
    bank_angle_plot(:,n) = bank_angle;
    
    n = n + 1;
end

% position plots
figure;
% plot 1
subplot(3,1,1)
plot(t_plot,p_plot(1,:))
title('north vs t')
xlabel 'time(sec)';
ylabel 'pos north (m)';

% plot 2
subplot(3,1,2)
plot(t_plot,p_plot(2,:))
title('east vs t')
xlabel 'time(sec)';
ylabel 'pos east (m)';

% plot 3
subplot(3,1,3)
plot(t_plot,p_plot(3,:))
title('down vs t')
xlabel 'time(sec)';
ylabel 'pos down (m)';



figure
subplot(2,1,1)
plot(t_plot, chi_plot*180/pi)
title('Heading angle vs t')
xlabel 'time(sec)';
ylabel 'chi (deg)';

subplot(2,1,2)
plot(t_plot, gamma_plot*180/pi)
title('Climb angle vs t')
xlabel 'time(sec)';
ylabel 'gamma (deg)';

figure
subplot(3,1,1)
plot(t_plot, L_plot)
title('Lift vs t')
xlabel 'time(sec)';
ylabel 'Lift (N)';

subplot(3,1,2)
plot(t_plot, D_plot)
title('Drag vs t')
xlabel 'time(sec)';
ylabel 'Drag (N)';

subplot(3,1,3)
plot(t_plot, T_plot)
title('Thrust vs t')
xlabel 'time(sec)';
ylabel 'Thrust (N)';


figure
subplot(2,1,1)
plot(t_plot, alpha_plot*180/pi)
title('angle of attack vs t')
xlabel 'time(sec)';
ylabel 'alpha (deg)';

subplot(2,1,2)
plot(t_plot, bank_angle_plot*180/pi)
title('bank angle vs t')
xlabel 'time(sec)';
ylabel 'bank angle (deg)';

figure;
plot3(p_plot(2,:),p_plot(1,:), -p_plot(3,:))
title('x vs y vs z')
xlabel 'east';
ylabel 'north';
zlabel 'neg_down';