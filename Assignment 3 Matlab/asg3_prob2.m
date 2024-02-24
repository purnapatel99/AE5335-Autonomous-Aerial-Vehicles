clear
clc

load asg3_prob2_coefficients.mat

m = 0.5;
L = 0.175;
J = diag([2.32 2.32 4]*1E-3);
Km = 1.5e-9;
Kf = 6.11e-8;


g = 9.81;

t_plot = [];
p_plot = [];
p_dot_plot = [];
attitude_plot = [];
w_plot = [];
F_plot = [];
M_b_plot = [];
motor_speed_plot = [];
n = 1;
q_plot = [];
for t = 0 : 0.1 : 60

    p = coeff_p(:,5) + coeff_p(:,4)*t  + coeff_p(:,3)*t^2 + coeff_p(:,2)*t^3 + coeff_p(:,1)*t^4;
    p_dot = coeff_p(:,4)  + 2*coeff_p(:,3)*t + 3*coeff_p(:,2)*t^2 + 4*coeff_p(:,1)*t^3;
    p_2dot = 2*coeff_p(:,3) + 6*coeff_p(:,2)*t + 12*coeff_p(:,1)*t^2;
    p_3dot = 6*coeff_p(:,2) + 24*coeff_p(:,1)*t;
    p_4dot = 24*coeff_p(:,1);


    psi = atan2(p_dot(2), p_dot(1));
    psi_dot = (p_dot(1)*p_2dot(2) - p_dot(2)*p_2dot(1))/(p_dot(1)^2 + p_dot(2)^2);
    psi_2dot = ((p_dot(1)^2 + p_dot(2)^2)*(p_dot(1)*p_3dot(2) - p_dot(2)*p_3dot(1))...
        - (p_dot(1)*p_2dot(2) - p_dot(2)*p_2dot(1))*(2*p_dot(1)*p_2dot(1) + 2*p_dot(2)*p_2dot(2)))...
        /(p_dot(1)^2 + p_dot(2)^2)^2;


    F_t = m*p_2dot - [0;0;m*g];
    
    F = norm(F_t);

    k_tb = -F_t/F;

    j_ta1 = [-sin(psi); cos(psi); 0];

    cross_jk = cross(j_ta1, k_tb);
    i_tb = cross_jk/norm(cross_jk);

    j_tb = cross(k_tb, i_tb);

    R_tb = [i_tb j_tb k_tb];

    temp_R = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1] * R_tb;

    th = atan2(-temp_R(3,1), temp_R(1,1));
    phi = atan2(-temp_R(2,3),temp_R(2,2));

    temp1 = (m*R_tb'*p_3dot);

    p_w = (1/F)*dot(temp1, [0;1;0]);
    q_w = (1/F)*dot(-temp1, [1;0;0]);

    temp2 = pinv([k_tb j_ta1 i_tb])*(R_tb*[p_w;q_w;0] - psi_dot*[0;0;1]);

    r_w = -temp2(1);
    th_dot = temp2(2);
    phi_dot = temp2(3);

    w_btb = [p_w;q_w;r_w];

    w_ttb = psi_dot*[0;0;1] + th_dot*j_ta1 + phi_dot*i_tb;
    
    F_dot = m^2*(p_2dot(1)*p_3dot(1) + p_2dot(2)*p_3dot(2) + (p_2dot(3)-g)*p_3dot(3)) / F;

    F_2dot = (m^2*(p_3dot(1)^2 + p_2dot(1)*p_4dot(1) + p_3dot(2)^2 +...
        p_2dot(2)*p_4dot(2) + p_3dot(3)^2 + (p_2dot(3)-g)*p_4dot(3)) - F_dot^2) / F;

    F_b = [0;0;-F];
    F_b_dot = [0;0;-F_dot];
    F_b_2dot = [0;0;-F_2dot];


    temp3 = (1/F)*(m*R_tb'*p_4dot - F_b_2dot - 2*cross(w_btb, F_b_dot) - cross(w_btb, cross(w_btb, F_b)));

    p_w_dot = temp3(2);
    q_w_dot = -temp3(1);

    i_tb_dot = [-psi_dot*sin(psi)*cos(th)-th_dot*cos(psi)*sin(th);
                -th_dot*sin(th)*sin(psi) + psi_dot*cos(th)*cos(psi);
                -th_dot*cos(th)];

    j_ta1_dot = psi_dot*[-cos(psi); -sin(psi); 0];

    temp4 = pinv([-k_tb j_ta1 i_tb])*(R_tb*[p_w_dot;q_w_dot;0] - psi_2dot*[0;0;1] - th_dot*j_ta1_dot - phi_dot*i_tb_dot);

    r_w_dot = temp4(1);

    w_btb_dot = [p_w_dot; q_w_dot; r_w_dot];

    M_b = J*w_btb_dot + cross(w_btb, J*w_btb);

    const_mat = [-Kf -Kf -Kf -Kf;
                 0 -Kf*L 0 Kf*L;
                 Kf*L 0 -Kf*L 0;
                 Km -Km Km -Km];
    motor_speed = sqrt(pinv(const_mat)*[-F;M_b]);

    t_plot(:,n) = t;
    attitude_plot(:,n) = [psi; th; phi];
    w_plot(:,n) = w_btb;
    p_plot(:,n) = p;
    p_dot_plot(:,n) = p_dot;
    F_plot(:,n) = F;
    M_b_plot(:,n) = M_b;
    motor_speed_plot(:,n) = motor_speed;

    qe = euler2quat([psi; th; phi]);
    q_plot(:,n) = qe;

    

    n = n + 1;

end

x_dynamics = [p_plot' p_dot_plot' attitude_plot' w_plot' q_plot'];
t_dynamics  = t_plot';
figure; %attitude angles plots
% plot 1
subplot(3,1,1)
plot(t_plot,attitude_plot(1,:)*180/pi)
title('Yaw vs t')
xlabel 'time(sec)';
ylabel 'psi (deg)';

% plot 2
subplot(3,1,2)
plot(t_plot,attitude_plot(2,:)*180/pi)
title('Pitch vs t')
xlabel 'time(sec)';
ylabel 'theta (deg)';

% plot 3
subplot(3,1,3)
plot(t_plot,attitude_plot(3,:)*180/pi)
title('Roll vs t')
xlabel 'time(sec)';
ylabel 'phi (deg)';


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

%velocity plot
figure;
% plot 1
subplot(3,1,1)
plot(t_plot,p_dot_plot(1,:))
title('vel north vs t')
xlabel 'time(sec)';
ylabel 'vel north (m/s)';

% plot 2
subplot(3,1,2)
plot(t_plot,p_dot_plot(2,:))
title('vel east vs t')
xlabel 'time(sec)';
ylabel 'vel east (m/s)';

% plot 3
subplot(3,1,3)
plot(t_plot,p_dot_plot(3,:))
title('vel down vs t')
xlabel 'time(sec)';
ylabel 'vel down (m/s)';


% %3D trajectory plot for visualition
% figure;
% plot3(p_plot(2,:),p_plot(1,:), -p_plot(3,:))
% title('phi vs t')
% xlabel 'east';
% ylabel 'north';
% zlabel 'neg_down';

% Omega plots
figure;
% plot 1
subplot(3,1,1)
plot(t_plot,w_plot(1,:)*180/pi)
title('omega_b p vs t')
xlabel 'time(sec)';
ylabel 'p (deg/s)';

% plot 2
subplot(3,1,2)
plot(t_plot,w_plot(2,:)*180/pi)
title('omega_b q vs t')
xlabel 'time(sec)';
ylabel 'q (deg/s)';

% plot 3
subplot(3,1,3)
plot(t_plot,w_plot(3,:)*180/pi)
title('omega_b r vs t')
xlabel 'time(sec)';
ylabel 'r (deg/s)';


% force plot
figure;
plot(t_plot,F_plot(1,:))
title('Total thrust vs t')
xlabel 'time(sec)';
ylabel 'F (N)';

figure
% plot 1
subplot(3,1,1)
plot(t_plot,M_b_plot(1,:))
title('Moment x vs t')
xlabel 'time(sec)';
ylabel 'Moment in x (Nm)';

% plot 2
subplot(3,1,2)
plot(t_plot,M_b_plot(2,:))
title('Moment y vs t')
xlabel 'time(sec)';
ylabel 'Moment in y (Nm)';

% plot 3
subplot(3,1,3)
plot(t_plot,M_b_plot(3,:))
title('Moment z vs t')
xlabel 'time(sec)';
ylabel 'Moment in z (Nm)';




% Motor speed plots
figure;
% plot 1
subplot(4,1,1)
plot(t_plot,motor_speed_plot(1,:))
title('motor 1 rotation speed vs t')
xlabel 'time(sec)';
ylabel 'speed motor 1';

% plot 2
subplot(4,1,2)
plot(t_plot,motor_speed_plot(2,:))
title('motor 2 rotation speed vs t')
xlabel 'time(sec)';
ylabel 'speed motor 2';

% plot 3
subplot(4,1,3)
plot(t_plot,motor_speed_plot(3,:))
title('motor 3 rotation speed vs t')
xlabel 'time(sec)';
ylabel 'speed motor 3';

% plot 4
subplot(4,1,4)
plot(t_plot,motor_speed_plot(4,:))
title('motor 4 rotation speed vs t')
xlabel 'time(sec)';
ylabel 'speed motor 4';