function motor_speeds = my_dynamics(t)

load asg3_prob2_coefficients.mat

m = 0.5;
L = 0.175;
J = diag([2.32 2.32 4]*1E-3);
Km = 1.5E-9;
Kf = 6.11E-8;

g = 9.81;

p = coeff_p(:,5) + coeff_p(:,4)*t  + coeff_p(:,3)*t^2 + coeff_p(:,2)*t^3 + coeff_p(:,1)*t^4;
p_dot = coeff_p(:,4)  + 2*coeff_p(:,3)*t + 3*coeff_p(:,2)*t^2 + 4*coeff_p(:,1)*t^3;
p_2dot = 2*coeff_p(:,3) + 6*coeff_p(:,2)*t + 12*coeff_p(:,1)*t^2;
p_3dot = 6*coeff_p(:,2) + 24*coeff_p(:,1)*t;
p_4dot = 24*coeff_p(:,1);

%     p = [0;0;-10];
%     p_dot = [0;0;0];
%     p_2dot = [0;0;0];
%     p_3dot = [0;0;0];
%     p_4dot = [0;0;0];

psi = atan2(p_dot(2), p_dot(1));
%     psi*180/pi
psi_dot = (p_dot(1)*p_2dot(2) - p_dot(2)*p_2dot(1))/(p_dot(1)^2 + p_dot(2)^2);
psi_2dot = ((p_dot(1)^2 + p_dot(2)^2)*(p_dot(1)*p_3dot(2) - p_dot(2)*p_3dot(1))...
    - (p_dot(1)*p_2dot(2) - p_dot(2)*p_2dot(1))*(2*p_dot(1)*p_2dot(1) + 2*p_dot(2)*p_2dot(2)))...
    /(p_dot(1)^2 + p_dot(2)^2)^2;

%     psi = 0;
%     psi_dot = 0;
%     psi_2dot = 0;

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
%     th = acos(temp_R(1,1));

phi = atan2(-temp_R(2,3),temp_R(2,2));
%     phi = acos(temp_R(2,2));

temp1 = (m*R_tb'*p_3dot);

p_w = (1/F)*dot(temp1, [0;1;0]);
q_w = (1/F)*dot(-temp1, [1;0;0]);

%   F_dot = dot(-temp1, [0;0;1])

temp2 = pinv([-k_tb j_ta1 i_tb])*(R_tb*[p_w;q_w;0] - psi_dot*[0;0;1]);

r_w = temp2(1);
th_dot = temp2(2);
phi_dot = temp2(3);

w_btb = [p_w;q_w;r_w];

w_ttb = psi_dot*[0;0;1] + th_dot*j_ta1 + phi_dot*i_tb;

F_dot = m^2*(p_2dot(1)*p_3dot(1) + p_2dot(2)*p_3dot(2) + (p_2dot(3)-g)*p_3dot(3)) / F;

F_2dot = (m^2*(p_3dot(1)^2 + p_2dot(1)*p_4dot(1) + p_3dot(2)^2 +...
    p_2dot(2)*p_4dot(2) + p_3dot(3)^2 + (p_2dot(3)-g)*p_4dot(3)) - F_dot^2) / F;

%     F_dot = 0;
%     F_2dot = 0;

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

% temp4 = pinv([-k_tb j_ta1 i_tb])*(R_tb*[p_w_dot;q_w_dot;0] - psi_2dot*[0;0;1]...
% - th_dot*(j_ta1_dot + cross(w_ttb,j_ta1)) - phi_dot*(i_tb_dot + cross(w_ttb, i_tb)));

r_w_dot = temp4(1);

w_btb_dot = [p_w_dot; q_w_dot; r_w_dot];

M_b = J*w_btb_dot + cross(w_btb, J*w_btb);

const_mat = [-Kf -Kf -Kf -Kf;
             0 -Kf*L 0 Kf*L;
             Kf*L 0 -Kf*L 0;
             Km -Km Km -Km];
motor_speeds = pinv(const_mat)*[-F;M_b];

end