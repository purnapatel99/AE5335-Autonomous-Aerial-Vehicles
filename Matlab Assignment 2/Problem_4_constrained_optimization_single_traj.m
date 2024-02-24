clear variables; close all; clc

%% Problem Data

p1 = [0 0]';
p2 = [10 -5]';
p3 = [30, 5]';
p4 = [40, 10]';
X1 = 10*pi/180;
X2 = 50*pi/180;
X3 = 120*pi/180;
X4 = 45*pi/180;

p = [p1, p2, p3];
X = [X1, X2, X3];

size_p = size(p);

%% Initial Guess
% Fitting a 3rd degree polynomil curve between the 2 points and stiching it.
% Varying x constantly with t and calculating the initial guesses for the coefficitnes. 
A_init = [1 p(1,1) p(1,1)^2 p(1,1)^3;
          1 p(1,2) p(1,2)^2 p(1,2)^3;
          1 p(1,3) p(1,3)^2 p(1,3)^3;
          0 1 2*p(1,1) 3*p(1,1)^2;
          0 1 2*p(1,2) 3*p(1,2)^2
          0 1 2*p(1,3) 3*p(1,3)^2];
b_init = [p(2,1) p(2,2) p(2,3) tan(X(1)) tan(X(2)) tan(X(3))]';

a = pinv(A_init)*b_init;
ax0 = p(1,1);
ax1 = p(1,3) - p(1,1);
ay0 = a(1) + a(2)*ax0 + a(3)*ax0^2 + a(4)*ax0^3;
ay1 = a(2)*ax1 + 2*a(3)*ax0*ax1 + 3*a(4)*ax0^2*ax1;
ay2 = a(3)*ax1^2 + 3*a(4)*ax0*ax1^2;
ay3 = a(4)*ax1^3;
t05 = (p(1,2) - p (1,1)) / (p(1,3) - p (1,1));
x_init = [ax0; ax1; 0; 0; 0; ay0; ay1; ay2; ay3; 0; t05];

% In non-trivial problems, a lot more thought needs to go into guess
% generation


%% Optimization

%----- Example of using 'optimoptions' to change the tool configuration
solver_options	= optimoptions('fmincon', 'MaxIterations', 1E4, ...
	'OptimalityTolerance', 1E-6, 'ConstraintTolerance', 1E-6, ...
	'MaxFunctionEvaluations', 1E6, 'Display', 'iter');

% Optimizing the trajectory using fmincon

[x_opt, fval, exit_flag, solver_output] = ....
    fmincon(@(x) my_cost(x, 100, 0, "acc"), ...
    x_init(:,1), [], [], [], [], [], [], ...									% <--- The two zeros here are lower bounds for the decision variables
    @(x) my_constraints(x, p(:,1), p(:,2), p(:,3), X(1), X(2), X(3)), ...
    solver_options);
% x_opt = [];
% for i = 1:(size_p(2)-1)
%     [x_opt(:,i), fval, exit_flag, solver_output] = ....
% 	    fmincon(@(x) my_cost(x, 100, 0, "acc"), ...
% 	    x_init(:,1), [], [], [], [], [], [], ...									% <--- The two zeros here are lower bounds for the decision variables
% 	    @(x) my_constraints(x, p(:,i), p(:,i+1), X(i), X(i+1)), ...
% 	    solver_options);
% 
% end
% 
% x_opt2 = [];
% for i = 1:(size_p(2)-1)
%     [x_opt2(:,i), fval, exit_flag, solver_output] = ....
% 	    fmincon(@(x) my_cost(x, 100, 1, "r_curvature"), ...
% 	    x_init(:,1), [], [], [], [], [], [], ...									% <--- The two zeros here are lower bounds for the decision variables
% 	    @(x) my_constraints(x, p(:,i), p(:,i+1), X(i), X(i+1)), ...
% 	    solver_options);
% 
% end

%% Generating datd and plotting the reults
arc_length(x_opt, 1000, 1)
l = @(t) arc_length(x_opt, 1000, t);
% kappa = @(t) curvature(x_opt, t);
xt = @(t) traj_x(x_opt,t);
yt = @(t) traj_y(x_opt,t);
% l2 = @(t) arc_length(x_opt2, 1000, t);
% kappa2 = @(t) curvature(x_opt2, t);
% xt2 = @(t) traj_x(x_opt2,t);
% yt2 = @(t) traj_y(x_opt2,t);
% 
figure
fplot(xt, yt, [0,1])
% xt2 = @(t) traj_x(x_opt2,t);
% yt2 = @(t) traj_y(x_opt2,t);
% 
% figure
% fplot(xt, yt, [0,2])
% hold on
% fplot(xt2, yt2, [0,2])
% 
% figure
% fplot(l, xt, [0,2])
% hold on
% fplot(l2, xt2, [0,2])
% 
% figure
% fplot(l, yt, [0,2])
% hold on
% fplot(l2, yt2, [0,2])
% 
% figure
% fplot(l, kappa, [0,2])
% hold on
% fplot(l2, kappa2, [0,2])









%% Function to calculate total arc length of the trajectory
function total_length = arc_length(x, M, t_end)
    dt = 1/M;
    length = 0;
 
    for i = 0 : M
        t = i*dt;
        if t >= t_end
            break;
        end
        dpx = x(2) + 2*x(3)*t + 3*x(4)*t^2 + 4*x(5)*t^3;
        dpy = x(7) + 2*x(8)*t + 3*x(9)*t^2 + 4*x(10)*t^3;
        d_length = sqrt(dpx^2 + dpy^2);

        if i == 0 || i == M
            length = length + 0.5*dt*d_length;
        else
            length = length + dt*d_length;
        end
    end
    total_length = length;
end 

%% These functions returns px and py of the trajectories of given coefficitnts
function px = traj_x(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    px = x(1,i+1) + x(2,i+1)*t + x(3,i+1)*t^2 + x(4,i+1)*t^3 + x(5,i+1)*t^4;
%     py = x(6,i+1) + x(7,i+1)*t + x(8,i+1)*t^2 + x(9,i+1)*t^3 + x(10,i+1)*t^4;
%     p = [px, py];
end

function py = traj_y(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
%     px = x(1,i+1) + x(2,i+1)*t + x(3,i+1)*t^2 + x(4,i+1)*t^3 + x(5,i+1)*t^4;
    py = x(6,i+1) + x(7,i+1)*t + x(8,i+1)*t^2 + x(9,i+1)*t^3 + x(10,i+1)*t^4;
%     p = [px, py];
end

%% Function to calculate curvature at given t
function kappa = curvature(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    dpx = x(2,i+1) + 2*x(3,i+1)*t + 3*x(4,i+1)*t^2 + 4*x(5,i+1)*t^3;
    dpy = x(7,i+1) + 2*x(8,i+1)*t + 3*x(9,i+1)*t^2 + 4*x(10,i+1)*t^3;
    ddpx = 2*x(3,i+1) + 6*x(4,i+1)*t + 12*x(5,i+1)*t^2;
    ddpy = 2*x(8,i+1) + 6*x(9,i+1)*t + 12*x(10,i+1)*t^2;
    dp = [dpx; dpy; 0];
    ddp = [ddpx; ddpy; 0];
    temp1 = norm(cross(dp, ddp));
    temp2 = norm(dp);
    kappa = temp1/temp2^3;
end

%% Cost function for optimization
function f_ = my_cost(x, M, alpha, method)
    if method == "r_curvature"
        alpha_cur = alpha;
        alpha_acc = 0;
    else
        alpha_cur = 0;
        alpha_acc = alpha;
    end

	dt = 1 / M;
    H_tilda = zeros(11,11);
    radius_curv = 0;
    for i = 0 : M
        t = i*dt;
        H_temp_1 = [0 1 2*t 3*t^2 4*t^3 0 0 0 0 0 0];
        H_temp_2 = [0 0 0 0 0 0 1 2*t 3*t^2 4*t^3 0];
        H_temp_3 = [0 0 2 6*t 12*t^2 0 0 0 0 0 0];
        H_temp_4 = [0 0 0 0 0 0 0 2 6*t 12*t^2 0];
        alpha1 = 1.0;
%         alpha = 1.0; % 0.0175
        H_i = alpha1*(H_temp_1'*H_temp_1 + H_temp_2'*H_temp_2) ...
            + alpha_acc*(H_temp_3'*H_temp_3 + H_temp_4'*H_temp_4);
        if i == 0 || i == M
            H_tilda = H_tilda + 0.5*dt*H_i;
%             radius_curv = radius_curv + 0.5*dt*(1/curvature(x, i*dt));
        else
            H_tilda = H_tilda + dt*H_i;
%             radius_curv = radius_curv + dt*(1/curvature(x, i*dt));
        end
    end

	f_ = x'*H_tilda*x + alpha_cur*radius_curv;	
end

%% Defining constraing for optimization
function [g_, h_] = my_constraints(x, p0, p1, p2, X0, X1, X2)
    M = 10;

    % Heading constrains at the end points of the trajectory
    A_g = [0 0 0 0 0 0 -sign(sin(X0)) 0 0 0 0;
           0 -sign(cos(X0)) 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 -sign(sin(X1)) -2*x(11)*sign(sin(X1)) -3*x(11)^2*sign(sin(X1)) -4*x(11)^3*sign(sin(X1)) 0;
           0 -sign(cos(X1)) -2*x(11)*sign(cos(X1)) -3*x(11)^2*sign(cos(X1)) -4*x(11)^3*sign(cos(X1)) 0 0 0 0 0 0;
           0 0 0 0 0 0 -sign(sin(X2)) -2*sign(sin(X2)) -3*sign(sin(X2)) -4*sign(sin(X2)) 0;
           0 -sign(cos(X2)) -2*sign(cos(X2)) -3*sign(cos(X2)) -4*sign(cos(X2)) 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 1;
           0 0 0 0 0 0 0 0 0 0 -1];
%            zeros(M-1,10)];

    % Veloity constrain for the trajectory
%     dt = 1/M;
%     vel = zeros(M+3,1);
%     for i = 1 : M-1
%         t = i*dt;
%         dpx = x(2) + 2*x(3)*t + 3*x(4)*t^2 + 4*x(5)*t^3;
%         dpy = x(7) + 2*x(8)*t + 3*x(9)*t^2 + 4*x(10)*t^3;
%         vel(i+4,1) = -sqrt(dpx^2 + dpy^2);
%     end

    
% 	g_ = A_g*x + vel + 1*[[sin(X0);cos(X0);sin(X1);cos(X1)];ones(M-1,1)];%A_g*x + 0*[1;1;1;1]; [sin(X0);cos(X0);sin(X1);cos(X1)]
	g_ = A_g*x + [1*[1;1;1;1;1;1];-0.9;0.1];
	% Defining equality constrain
    A = [1 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 0 0 0 0;
         1 x(11) x(11)^2 x(11)^3 x(11)^4 0 0 0 0 0 0;
         0 0 0 0 0 1 x(11) x(11)^2 x(11)^3 x(11)^4 0;
         1 1 1 1 1 0 0 0 0 0 0;
         0 0 0 0 0 1 1 1 1 1 0;
         0 tan(X0) 0 0 0 0 -1 0 0 0 0;
         0 tan(X1) 2*x(11)*tan(X1) 3*x(11)^2*tan(X1) 4*x(11)^3*tan(X1) 0 -1 -2*x(11) -3*x(11)^2 -4*x(11)^3 0
         0 tan(X2) 2*tan(X2) 3*tan(X2) 4*tan(X2) 0 -1 -2 -3 -4 0];
    b = [p0(1) p0(2) p1(1) p1(2) p2(1) p2(2) 0 0 0]';
	h_ = A*x - b;
	
end