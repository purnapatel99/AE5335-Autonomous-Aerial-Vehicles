%{
SOFTWARE LICENSE
----------------
Copyright (c) 2021 by Raghvendra V. Cowlagi

Permission is hereby granted to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in
the Software, including the rights to use, copy, modify, merge, copies of
the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:  

* The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.
* The Software, and its copies or modifications, may not be distributed,
published, or sold for profit. 
* The Software, and any substantial portion thereof, may not be copied or
modified for commercial or for-profit use.

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other
dealings in the software.      


PROGRAM DESCRIPTION
-------------------
This program illustrates the implementation of an extended Kalman filter
(EKF) for attitude estimation. Attitude is described by z-y-x aircraft
Euler angles. Measurements considered are of rate gyros and a fictitous
sensor that directly measures all Euler angles. Measurement bias is not
considered.
%}

clear variables; close all; clc;

load ae5335_ekfdemo_nobias_synthetic_data.mat 

n_pts = numel(timestamps);

R_me_covariance		= (std_dev_compass^2) * eye(3);							% Measurement error covariance (3x3 matrix)
P_ee_covariance     = R_me_covariance;										% Estimation error covariance
Q_pe_covariance		= (std_dev_rategyro^2) * eye(3);						% Process noise covariance

% Initial estimate same as direct measurement; zero initial estimate of all biases
x_hat				= z_compass3D_mb(:, 1);
x_int				= z_compass3D_mb(:, 1);

% Store mean and e.e. covariance at each time step
x_hat_store			= zeros(3, n_pts);
trace_P_store		= zeros(1, n_pts);
x_hat_store(:, 1)	= x_hat(1:3);
trace_P_store(1)	= trace(P_ee_covariance); 

% Store Euler angles calculated by integrating attitude kinematics (i.e., no filter)
x_int_store			= zeros(3, n_pts);
x_int_store(:, 1)	= x_int;

for k = 2:n_pts
    dt = timestamps(k) - timestamps(k - 1);
	
	% Physics-based predictive model    
	e_thta	= x_hat(2);
	e_phi	= x_hat(3);
	
	zk		= z_compass3D_mb(:, k);	% New magnetometer measurements
	
	% Predictive estimate (RK4) and e.e. covariance
	u_k1	= z_rategyro_rb(:, k - 1);
	u_k		= z_rategyro_rb(:, k);
	u12		= 0.5*(u_k1 + u_k);												% Approx (p,q,r) at (t + 0.5dt)
	k1		= dt*euler_kinematics(x_hat(1:3),				u_k1);
	k2		= dt*euler_kinematics((x_hat(1:3) + 0.5*k1),	u12);
	k3		= dt*euler_kinematics((x_hat(1:3) + 0.5*k2),	u12);
	k4		= dt*euler_kinematics((x_hat(1:3) + k3),		u_k);
	x_minus	= x_hat + ((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4);			% The predictive state update can use nonlinear kinematic equations

	
	% Linearization
	F0		= [...
		0	(u_k1(2)*sin(e_phi) + u_k1(3)*cos(e_phi))*tan(e_thta)*sec(e_thta)	(u_k1(2)*cos(e_phi) - u_k1(3)*sin(e_phi))*sec(e_thta); ...
		0	0																-(u_k1(2)*sin(e_phi) + u_k1(3)*cos(e_phi)); ...
		0	(u_k1(2)*sin(e_phi) + u_k1(3)*cos(e_phi))*sec(e_thta)^2				(u_k1(2)*cos(e_phi) - u_k1(3)*sin(e_phi))*tan(e_thta)];
	
    H0		= (1 / cos(e_thta) ) * [0 sin(e_phi) cos(e_phi); ...
		0 cos(e_phi)*cos(e_thta) -sin(e_phi)*cos(e_thta); ...
		cos(e_thta) sin(e_phi)*sin(e_thta) cos(e_phi)*sin(e_thta)];
	
	F_k1	= eye(3) + F0*dt;
	G2_k1	= -H0*dt;
	
	P_minus= F_k1*P_ee_covariance*F_k1' + G2_k1*Q_pe_covariance*G2_k1';		% The e.e. covariance update uses linearization...
	
	% Measurement model
    C	= eye(3);
	
	% Kalman gain
    L_Kalman_gain	= P_minus*C' / (C*P_minus*C' + R_me_covariance);		%... as does the Kalman gain computation
	
	% Recursive WLSE method to update predictive estimate with new
	% measurement
	x_hat	= x_minus + L_Kalman_gain * ( zk - C*x_minus );					% New estimate
	P_ee_covariance		= (eye(3) - L_Kalman_gain * C) * P_minus;			% New e.e. covariance, again linearized C 
	
	% Just for comparison: what would happen if we only predicted, did not "correct"
	k1	= dt*euler_kinematics(x_int(1:3),				u_k1);
	k2	= dt*euler_kinematics((x_int(1:3) + 0.5*k1),	u12);
	k3	= dt*euler_kinematics((x_int(1:3) + 0.5*k2),	u12);
	k4	= dt*euler_kinematics((x_int(1:3) + k3),		u_k);
	x_int	= x_int + ((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4); 
	
	% Record
	x_hat_store(:, k)	= x_hat(1:3);
	x_int_store(:, k)	= x_int(1:3);
	trace_P_store(k)	= trace(P_ee_covariance);
	% The trace of estimation error covariance should decrease with time, 
	% indicating increasing confidence in the state estimate
end

figsize = [0, 0.04, 0.8, 0.8];
figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);

subplot(311); plot(timestamps, euler_angles_true(1,:)*180/pi, 'LineWidth', 2); hold on;
subplot(311); plot(timestamps, x_hat_store(1,:)*180/pi, 'LineWidth', 2); 
subplot(311); plot(timestamps, x_int_store(1,:)*180/pi, '--', 'LineWidth', 2); 
subplot(311); plot(timestamps, z_compass3D_mb(1,:)*180/pi, 'LineWidth', 0.5); 
xlabel('Time (s)'); ylabel('\psi (deg)'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
legend('True', 'EKF estimate', 'Prediction only', 'Measurement only')
% legend('True', 'EKF estimate', 'Prediction only')

subplot(312); plot(timestamps, euler_angles_true(2,:)*180/pi, 'LineWidth', 2); hold on;
subplot(312); plot(timestamps, x_hat_store(2,:)*180/pi, 'LineWidth', 2);
subplot(312); plot(timestamps, x_int_store(2,:)*180/pi, '--', 'LineWidth', 2); 
subplot(312); plot(timestamps, z_compass3D_mb(2,:)*180/pi, 'LineWidth', 0.5); 
xlabel('Time (s)'); ylabel('\theta (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
legend('True', 'EKF estimate', 'Prediction only', 'Measurement only')
% legend('True', 'EKF estimate', 'Prediction only')

subplot(313); plot(timestamps, euler_angles_true(3,:)*180/pi, 'LineWidth', 2); hold on;
subplot(313); plot(timestamps, x_hat_store(3,:)*180/pi, 'LineWidth', 2);
subplot(313); plot(timestamps, x_int_store(3,:)*180/pi, '--', 'LineWidth', 2); 
subplot(313); plot(timestamps, z_compass3D_mb(3,:)*180/pi, 'LineWidth', 0.5); 
xlabel('Time (s)'); ylabel('\phi (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
legend('True', 'EKF estimate', 'Prediction only', 'Measurement only')
% legend('True', 'EKF estimate', 'Prediction only')

% exportgraphics(gcf, 'demo_ekf_attitude_nobias_2.png')


figsize = [0, 0.04, 0.8, 0.8];
figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);
plot(timestamps, trace_P_store, 'LineWidth', 4);
xlabel('Time (s)'); ylabel('trace(P)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

exportgraphics(gcf, 'demo_ekf_attitude_nobias_traceP.png')

function x_dot = euler_kinematics(x, pqr)
	e_thta_	= x(2);
	e_phi_	= x(3);

	x_dot	= [...
		-sin(e_thta_) 0 1; ...
		sin(e_phi_)*cos(e_thta_) cos(e_phi_) 0; ...
		cos(e_phi_)*cos(e_thta_) -sin(e_phi_) 0] \ pqr;
end