clear variables; close all; clc

%% Problem Data

p1 = [0 0]';
p2 = [10 -5]';
p3 = [30, 5]';

X1 = 10*pi/180;
X2 = 50*pi/180;
X3 = 120*pi/180;


p = [p1, p2, p3];
X = [X1, X2, X3];

size_p = size(p);

%% Initial Guess
% Fitting a 3rd degree polynomil curve between the 2 points and stiching it.
% Varying x constantly with t and calculating the initial guesses for the coefficitnes. 
x_init = [];
for i = 1:(size_p(2)-1)

    A_init = [1 p(1,i) p(1,i)^2 p(1,i)^3;
              1 p(1,i+1) p(1,i+1)^2 p(1,i+1)^3;
              0 1 2*p(1,i) 3*p(1,i)^2;
              0 1 2*p(1,i+1) 3*p(1,i+1)^2];
    b_init = [p(2,i) p(2,i+1) tan(X(i)) tan(X(i+1))]';
    
    a = pinv(A_init)*b_init;
    ax0 = p(1,i);
    ax1 = p(1,i+1) - p(1,i);
    ay0 = a(1) + a(2)*ax0 + a(3)*ax0^2 + a(4)*ax0^3;
    ay1 = a(2)*ax1 + 2*a(3)*ax0*ax1 + 3*a(4)*ax0^2*ax1;
    ay2 = a(3)*ax1^2 + 3*a(4)*ax0*ax1^2;
    ay3 = a(4)*ax1^3;
    x_init(:,i) = [ax0; ax1; 0; 0; 0; ay0; ay1; ay2; ay3; 0];
end


%% Optimization


solver_options	= optimoptions('fmincon', 'MaxIterations', 1E4, ...
	'OptimalityTolerance', 1E-6, 'ConstraintTolerance', 1E-6, ...
	'MaxFunctionEvaluations', 1E6, 'Display', 'iter');

% Optimizing the trajectory using fmincon
x_opt = [];
for i = 1:(size_p(2)-1)

    [x_opt(:,i), fval, exit_flag, solver_output] = ....
	    fmincon(@(x) my_cost(x, 100, 1), ...
	    x_init(:,1), [], [], [], [], [], [], ...									% <--- The two zeros here are lower bounds for the decision variables
	    @(x) my_constraints(x, p(:,i), p(:,i+1), X(i), X(i+1), 8), ...
	    solver_options);

end


%% Generating datd and plotting the reults
path_lenght1 = arc_length(x_opt, 1000, 2);
l = @(t) arc_length(x_opt, 1000, t);
kappa = @(t) curvature(x_opt, t);
xt = @(t) traj_x(x_opt,t);
yt = @(t) traj_y(x_opt,t);
angle = @(t) traj_angle(x_opt, t);
velocity = @(t) traj_vel(x_opt, t);


% figure
% fplot(xt, yt, [0,2])
% title('py vs px')
% xlabel 'px (m)';
% ylabel 'py (m)';
% hold on
% plot(p1(1), p1(2), '.', 'markersize', 20);
% plot(p2(1), p2(2), '.', 'markersize', 20);
% plot(p3(1), p3(2), '.', 'markersize', 20);
% 
% figure
% fplot(l, xt, [0,2])
% title('px vs path length')
% xlabel 'path length (m)';
% ylabel 'px (m)';
% 
% figure
% fplot(l, yt, [0,2])
% title('py vs path length')
% xlabel 'path length (m)';
% ylabel 'py (m)';
% 
% figure
% fplot(l, kappa, [0,2])
% title('curvature vs path length')
% xlabel 'path length (m)';
% ylabel 'curvature';
% 
% figure
% fplot(l, angle, [0,2])
% title('angle vs path length')
% xlabel 'path length (m)';
% ylabel 'angle (deg)';








%% Function to calculate total arc length of the trajectory
function total_length = arc_length(x, M, t_end)
    dt = 1/M;
    length = 0;
    
    size_x = size(x);
    for i = 1 : size_x(2)
        for j = 0 : M
            t = j*dt;
            if t+i-1 >= t_end
                break;
            end
            dpx = x(2,i) + 2*x(3,i)*t + 3*x(4,i)*t^2 + 4*x(5,i)*t^3;
            dpy = x(7,i) + 2*x(8,i)*t + 3*x(9,i)*t^2 + 4*x(10,i)*t^3;
            d_length = sqrt(dpx^2 + dpy^2);
    
            if j == 0 || j == M
                length = length + 0.5*dt*d_length;
            else
                length = length + dt*d_length;
            end
            
        end
        if t+i-1 >= t_end
            break;
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
end

function py = traj_y(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    py = x(6,i+1) + x(7,i+1)*t + x(8,i+1)*t^2 + x(9,i+1)*t^3 + x(10,i+1)*t^4;
end

function X = traj_angle(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    dpx = x(2,i+1) + 2*x(3,i+1)*t + 3*x(4,i+1)*t^2 + 4*x(5,i+1)*t^3;
    dpy = x(7,i+1) + 2*x(8,i+1)*t + 3*x(9,i+1)*t^2 + 4*x(10,i+1)*t^3;
    X = atan2(dpy, dpx)*180/pi;
end

function X = traj_vel(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    dpx = x(2,i+1) + 2*x(3,i+1)*t + 3*x(4,i+1)*t^2 + 4*x(5,i+1)*t^3;
    dpy = x(7,i+1) + 2*x(8,i+1)*t + 3*x(9,i+1)*t^2 + 4*x(10,i+1)*t^3;
    X = sqrt(dpy^2 + dpx^2);
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
function f_ = my_cost(x, M, u)

	dt = 1 / M;
    H_tilda = zeros(10,10);
    for i = 0 : M
        t = i*dt;
        H_temp_1 = [0 1 2*t 3*t^2 4*t^3 0 0 0 0 0];
        H_temp_2 = [0 0 0 0 0 0 1 2*t 3*t^2 4*t^3];
        H_temp_3 = [0 0 2 6*t 12*t^2 0 0 0 0 0];
        H_temp_4 = [0 0 0 0 0 0 0 2 6*t 12*t^2];
        
        H_i = (H_temp_1'*H_temp_1 + H_temp_2'*H_temp_2) ...
            + u*(H_temp_3'*H_temp_3 + H_temp_4'*H_temp_4);
        if i == 0 || i == M
            H_tilda = H_tilda + 0.5*dt*H_i;
        else
            H_tilda = H_tilda + dt*H_i;
        end
    end

	f_ = x'*H_tilda*x;	
end

%% Defining constraing for optimization
function [g_, h_] = my_constraints(x, p0, p1, X0, X1, v)
    M = 1000;

    % Heading constrains at the end points of the trajectory
    A_g = [0 0 0 0 0 0 -sign(sin(X0)) 0 0 0;
           0 -sign(cos(X0)) 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 -sign(sin(X1)) -2*sign(sin(X1)) -3*sign(sin(X1)) -4*sign(sin(X1));
           0 -sign(cos(X1)) -2*sign(cos(X1)) -3*sign(cos(X1)) -4*sign(cos(X1)) 0 0 0 0 0;
           zeros(M-1,10)];

    % Veloity constrain for the trajectory
    dt = 1/M;
    vel = zeros(M+3,1);
    for i = 1 : M-1
        t = i*dt;
        dpx = x(2) + 2*x(3)*t + 3*x(4)*t^2 + 4*x(5)*t^3;
        dpy = x(7) + 2*x(8)*t + 3*x(9)*t^2 + 4*x(10)*t^3;
        vel(i+4,1) = -sqrt(dpx^2 + dpy^2);
    end

    
	g_ = A_g*x + vel + [zeros(4,1);v*ones(M-1,1)];

    % Equality constrains
    A = [1 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 0 0 0;
         1 1 1 1 1 0 0 0 0 0;
         0 0 0 0 0 1 1 1 1 1;
         0 tan(X0) 0 0 0 0 -1 0 0 0;
         0 tan(X1) 2*tan(X1) 3*tan(X1) 4*tan(X1) 0 -1 -2 -3 -4];
    b = [p0(1) p0(2) p1(1) p1(2) 0 0]';
	h_ = A*x - b;
	
end