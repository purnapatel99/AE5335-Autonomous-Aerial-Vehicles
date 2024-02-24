%{
SOFTWARE LICENSE
----------------
Copyright (c) 2022 by Raghvendra V. Cowlagi

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
This program is a low-altitude quadrotor simulator that solves rigid body
equations with quadrotor forces and moments. The rotor spin rates can be
provided as inputs (see instructions below). Wind is not considered.
%}

%%
clear variables; close all; clc;
asg3_prob2
%% Simulation settings and aircraft parameters

t_final			= 60;

a0				= 0;
Earth.g			= 9.81;						% m/s^2
Earth.rho		= 1.225;					% kg/m^3

aircraft.Mass	= 0.5;						% kg
aircraft.MOI	= diag([2.32 2.32 4]*1E-3); % kg/m^2
aircraft.L		= 0.175;					% m
aircraft.k_F	= 6.11E-8;					% N/rpm^2
aircraft.k_M	= 1.5E-9;					% Nm/rpm^2
aircraft.Weight	= aircraft.Mass * Earth.g;

%----- Feedforward controller examples
% my_controller.ff= @ffcontrol_hover;				% Feedforward hover controller
% my_controller.ff= @ffcontrol_north;					% Feedforward hover controller
my_controller.ff= @ffcontrol_north_lookup;			% Feedforward hover controller (lookup table implementation)
%{
 NOTE: CONTRLLER FUNCTIONS RETURN THE SQUARES OF SPIN RATES.

 WRITE YOUR OWN BASED ON THESE EXAMPLES.
 ffcontrol_north_table PROVIDES AN EXAMPLE OF GETTING CONTROL INPUTS
 FROM A LOOKUP TABLE. DEFINE A LOOKUP TABLE WITH FIRST ROW OF TIMESTAMPS
 AND SECOND-FIFTH ROWS OF CONTROL VALUES (SQUARES OF SPIN RATES).

 DEVELOP A LOOKUP TABLE LIKE THIS FOR A GIVEN DESIRED TRAJECTORY USING THE
 DIFFERENTIAL FLATNESS PROPERTY. 
%}

my_control_lookup(1, :)		= linspace(0, t_final, 1001);					% Time stamps
for n = 1:1001
	t = my_control_lookup(1, n);
	my_control_lookup(2:5, n) = my_dynamics(t);
end
% my_control_lookup(2:5, :)	= ( ...											% Constant spin rates for this example
% 	(0.25)*sqrt( aircraft.Weight^2 + (aircraft.Mass*a0)^2 ) / ...
% 	aircraft.k_F )*ones(4, 100);
aircraft.control_lookup		= my_control_lookup;

%----- Feedback controller (none), do not change
my_controller.fb= @fbcontrol_zero;			% Zero feedback controller (placeholder)

%% Initial conditions
% CHANGE THESE INITIAL CONDITIONS AS NEEDED



p_I_init= [-32.1137126785674; -21.2816075411487; 0.677518190021755];		% Earth-relative position in inertial coordinates
Vg_I	= [12.4510; 13.9500; -0.1571];						                % Earth-relative velocity in body-fixed coordinates
e_psi	= 0.842114549614834; %0;											% Euler angle psi (yaw)
e_theta	= 0.192979010571882; %-atan(a0 / Earth.g);							% Euler angle theta (pitch)
e_phi	= -0.0960421388972161; %0;											% Euler angle phi (roll)
omega_B	= [0.00177442806699433; -0.00864912276424138; -0.0519087474532719];	% Inertial angular velocity in body-fixed coordinates
qe		= euler2quat([e_psi; e_theta; e_phi]);

x_init	= [p_I_init; Vg_I; e_psi; e_theta; e_phi; omega_B; qe];

%% Run simulation

[t_sim, x_sim] = ode45(@(t, x) ...
	quadrotor_dynamics(t, x, aircraft, my_controller, []), ...
	[0 t_final], x_init);
ffc_sim = zeros( length(t_sim), 4);
for m1 = 1:length(t_sim)
	u_ = my_controller.ff(t_sim(m1), aircraft);
	ffc_sim(m1, :) = u_';
end

%% Plot results
my_titles = {...
	'$pos_x $vs $time$'; '$pos_y $vs $time$'; '$pos_z $vs $time$'; ...
	'$vel_x $vs $time$'; '$vel_y $vs $time$'; '$vel_z $vs $time$'; ...
	'$yaw $vs $time$'; '$pitch $vs $time$'; '$roll $vs $time$'; ...
	'$omega_p $vs $time$'; '$omega_q $vs $time$'; '$omega_r $vs $time$'; ...
	'$q_0 $vs $time$'; '$q_1 $vs $time$'; '$q_2 $vs $time$'; '$q_3 $vs $time$'};
my_labels = {...
	'$p_x$ (m)'; '$p_y$ (m)'; '$p_z$ (m)'; ...
	'$V_x$ (m/s)'; '$V_y$ (m/s)'; '$V_z$ (m/s)'; ...
	'$\psi$ (deg)'; '$\theta$ (deg)'; '$\phi$ (deg)'; ...
	'$p$ (deg/s)'; '$q$ (deg/s)'; '$r$ (deg/s)'; ...
	'$q_0$'; '$q_1$'; '$q_2$'; '$q_3$'};
my_plot_scale = [ones(6,1); 180/pi*ones(6,1); ones(4,1)];
figure('Units', 'normalized', 'Position', [0.05 0.45 0.6 0.5], 'Name', 'States');
for m1 = 1:3 % 16 states incl. quaternions and Euler angles (redundant)
	subplot(3,1,m1);
	plot(t_sim, x_sim(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    hold on
    plot(t_dynamics, x_dynamics(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    title(my_titles{m1}, 'Interpreter', 'latex');
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end
figure('Units', 'normalized', 'Position', [0.05 0.45 0.6 0.5], 'Name', 'States');
for m1 = 4:6 % 16 states incl. quaternions and Euler angles (redundant)
	subplot(3,1,m1-3);
	plot(t_sim, x_sim(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    hold on
    plot(t_dynamics, x_dynamics(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    title(my_titles{m1}, 'Interpreter', 'latex');
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end
figure('Units', 'normalized', 'Position', [0.05 0.45 0.6 0.5], 'Name', 'States');
for m1 = 7:9 % 16 states incl. quaternions and Euler angles (redundant)
	subplot(3,1,m1-6);
	plot(t_sim, x_sim(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    hold on
    plot(t_dynamics, x_dynamics(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    title(my_titles{m1}, 'Interpreter', 'latex');
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end
figure('Units', 'normalized', 'Position', [0.05 0.45 0.6 0.5], 'Name', 'States');
for m1 = 10:12 % 16 states incl. quaternions and Euler angles (redundant)
	subplot(3,1,m1-9);
	plot(t_sim, x_sim(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    hold on
    plot(t_dynamics, x_dynamics(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    title(my_titles{m1}, 'Interpreter', 'latex');
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end
figure('Units', 'normalized', 'Position', [0.05 0.45 0.6 0.5], 'Name', 'States');
for m1 = 13:16 % 16 states incl. quaternions and Euler angles (redundant)
	subplot(2,2,m1-12);
	plot(t_sim, x_sim(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    hold on
    plot(t_dynamics, x_dynamics(:, m1)*my_plot_scale(m1), 'LineWidth', 1);
    title(my_titles{m1}, 'Interpreter', 'latex');
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end

my_labels = {...
	'$\Omega_1$ (rpm)'; '$\Omega_2$ (rpm)'; ...
	'$\Omega_3$ (rpm)'; '$\Omega_4$ (rpm)'};
figure('Units', 'normalized', 'Position', [0.65 0.05 0.3 0.5], 'Name', 'Controls');
for m1 = 1:4
	subplot(4,1,m1);
	plot(t_sim, ffc_sim(:, m1).^0.5, 'LineWidth', 2);
	ylabel(my_labels{m1}, 'Interpreter', 'latex'); xlabel('Time (s)');
	ax = gca;
	ax.FontName = 'Times New Roman';
	ax.FontSize = 14;
end

%% Quadrotor dynamics

function x_dot = quadrotor_dynamics(t_, x_, aircraft_, controller_, Earth_) 

p_I		= x_(1:3);															% Earth-relative position in inertial coordinates
Vg_I	= x_(4:6);															% Earth-relative velocity in body-fixed coordinates
psi_	= x_(7);															% Euler angle psi (yaw)
theta_	= x_(8);															% Euler angle theta (pitch)
phi_	= x_(9);															% Euler angle phi (roll)
omega_B	= x_(10:12);														% Inertial angular velocity in body-fixed coordinates
qe_		= x_(13:16);														% Quaternions

uff_	= controller_.ff(t_, aircraft_);									% Feedforward controller
ufb_	= controller_.fb(t_, x_, aircraft_);								% Feedforward controller
u_squared_spin_rates = uff_ + ufb_;		% rpm^2 units

%----- This is the Euler angle representation of the I-> B rotation matrix
% Body-fixed axis system convention is the same as that for fixed-wing
% aircraft, i.e., z points down.
Rot_IB	= ...
	[1 0 0; 0 cos(phi_) sin(phi_); 0 -sin(phi_) cos(phi_)] * ...
	[cos(theta_) 0 -sin(theta_); 0 1 0; sin(theta_) 0 cos(theta_)] * ...
	[cos(psi_) sin(psi_) 0; -sin(psi_) cos(psi_) 0; 0 0 1];
Rot_BI	= Rot_IB';

H_321	= [...
	-sin(theta_)			0			1; ...
	sin(phi_)*cos(theta_)	cos(phi_)	0; ...
	cos(phi_)*cos(theta_)	-sin(phi_)	0];

p_		= omega_B(1);
q_		= omega_B(2);
r_		= omega_B(3);

%----- This model can work with either Euler angles or Euler-Rodrigues
% parameters (quaternions). Both are included in the state variable, but
% only one needs to be used for I -> B transformations. Comment out of one
% these. 
%----- Translational kinematics and weight (using Euler angles)
dpos_dt	= Vg_I;

%----- Attitude kinematics (Euler angles)
deu_dt	= H_321 \ omega_B; 

%----- Attitude kinematics (Quaternions)
dq_dt	= 0.5*[0 -p_ -q_ -r_; p_ 0 r_ -q_; q_ -r_ 0 p_; r_ q_ -p_ 0]*qe_;

%----- Forces and moments
F_thrust_B	= [0; 0; -aircraft_.k_F*ones(1, 4)*u_squared_spin_rates];		% Rotor thrust in body-fixed coordinates
F_weight_I	= [0; 0; aircraft_.Weight];										% Weight in inertial coordinates
M_B		= [...
	0,	-aircraft_.k_F*aircraft_.L,	0,	aircraft_.k_F*aircraft_.L; ...
	aircraft_.k_F*aircraft_.L,	0,	-aircraft_.k_F*aircraft_.L,	0; ...
	aircraft_.k_M,	-aircraft_.k_M,	aircraft_.k_M,	-aircraft_.k_M] * ...	% External moments about CM in body-fixed coordinates
	u_squared_spin_rates;

%----- Translational dynamics (in inertial coordinates, using Euler angles)
dVg_dt	= (Rot_BI*F_thrust_B + F_weight_I) / aircraft_.Mass;

%----- Attitude dynamics
do_dt	= aircraft_.MOI \ ( M_B - cross(omega_B, aircraft_.MOI*omega_B) );

%----- Concatenate all state derivatives
x_dot	= [dpos_dt; dVg_dt; deu_dt; do_dt; dq_dt];

end

%% Controllers

function u_ = ffcontrol_hover(t_, aircraft_)
	u_ = (0.25)*aircraft_.Weight / aircraft_.k_F;
	% This is u_ = (\Omega_1^2, ... , \Omega_4^2)
end

function u_ = ffcontrol_north(t_, aircraft_)
	a0 = 1;
	u_ = (0.25)*sqrt( aircraft_.Weight^2 + (aircraft_.Mass*a0)^2 ) / aircraft_.k_F;
	% This is u_ = (\Omega_1^2, ... , \Omega_4^2)
end

function u_ = ffcontrol_north_lookup(t_, aircraft_)
	lookup_time_= aircraft_.control_lookup(1, :);
	t_indx		= find(lookup_time_ >= t_, 1, 'first');	
	if isempty(t_indx), t_indx = length(lookup_time_); end

	u_	= aircraft_.control_lookup(2:5, t_indx);
	% This is u_ = (\Omega_1^2, ... , \Omega_4^2)
end


function u_ = fbcontrol_zero(t_, x_, aircraft_)
	u_ = zeros(4, 1);
end