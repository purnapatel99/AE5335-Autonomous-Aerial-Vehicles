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
This program illustrates the use of the 'fmincon' tool for constrained
optimization. The problem considered is of maximizing the volume of a
cylinder subject to constraints on its height and surface area.
%}

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
    x_init(:,i) = [ax0; ax1; 0; 0; ay0; ay1; ay2; ay3];
end

% In non-trivial problems, a lot more thought needs to go into guess
% generation


%% Optimization

%----- Example of using 'optimoptions' to change the tool configuration
solver_options	= optimoptions('fmincon', 'MaxIterations', 1E4, ...
	'OptimalityTolerance', 1E-6, 'ConstraintTolerance', 1E-6, ...
	'MaxFunctionEvaluations', 1E6, 'Display', 'iter');

%----- Application of 'fmincon'
% Linear constraints and bounds can be encoded in the fmincon parameters.
% See the documentation. A separate function needs to be written for
% nonlinear constraints.
x_opt = [];
for i = 1:(size_p(2)-1)
    [x_opt(:,i), fval, exit_flag, solver_output] = ....
	    fmincon(@(x) my_cost(x, 1000, []), ...
	    x_init(:,1), [], [], [], [], [0; 0], [], ...									% <--- The two zeros here are lower bounds for the decision variables
	    @(x) my_constraints(x, p(:,i), p(:,i+1), X(i), X(i+1)), ...
	    solver_options);

%     fprintf('\n----- Optimal solution ----- \n')
%     fprintf('\t Optimal path coeff\t = %6.3f \n ', x_opt(:,i))
end
% arc_length(x_opt, 1000, 2)
curvature(x_opt, 2)
l = @(t) arc_length(x_opt, 1000, t);
kappa = @(t) curvature(x_opt, t);
% traj_x(x_opt, 2)
xt = @(t) traj_x(x_opt,t);
yt = @(t) traj_y(x_opt,t);
figure
fplot(xt, yt, [0,2])
figure
fplot(l, xt, [0,2])
figure
fplot(l, yt, [0,2])
figure
fplot(l, kappa, [0,1.99])










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
            dpx = x(2,i) + 2*x(3,i)*t + 3*x(4,i)*t^2;
            dpy = x(6,i) + 2*x(7,i)*t + 3*x(8,i)*t^2;
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


function px = traj_x(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    px = x(1,i+1) + x(2,i+1)*t + x(3,i+1)*t^2 + x(4,i+1)*t^3;
    py = x(5,i+1) + x(6,i+1)*t + x(7,i+1)*t^2 + x(8,i+1)*t^3;
    p = [px, py];
end

function py = traj_y(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    px = x(1,i+1) + x(2,i+1)*t + x(3,i+1)*t^2 + x(4,i+1)*t^3;
    py = x(5,i+1) + x(6,i+1)*t + x(7,i+1)*t^2 + x(8,i+1)*t^3;
    p = [px, py];
end


function kappa = curvature(x, t)
    size_x = size(x);
    t=min(max(t,0),size_x(2));
    i = t - mod(t,1);
    i = min(i, size_x(2) - 1);
    t = t-i;
    dpx = x(2,i+1) + 2*x(3,i+1)*t + 3*x(4,i+1)*t^2;
    dpy = x(6,i+1) + 2*x(7,i+1)*t + 3*x(8,i+1)*t^2;
    ddpx = 2*x(3,i+1) + 6*x(4,i+1)*t;
    ddpy = 2*x(7,i+1) + 6*x(8,i+1)*t;
    dp = [dpx; dpy; 0];
    ddp = [ddpx; ddpy; 0];
    temp1 = norm(cross(dp, ddp));
    temp2 = norm(dp);
    kappa = temp1/temp2^3;
end

function f_ = my_cost(x, M, parameter_2)
	dt = 1 / M;
    H_tilda = zeros(8,8);

    for i = 0 : M
        H_temp_1 = [0 1 2*i*dt 3*(i*dt)^2 0 0 0 0];
        H_temp_2 = [0 0 0 0 0 1 2*i*dt 3*(i*dt)^2];
        H_i = H_temp_1'*H_temp_1 + H_temp_2'*H_temp_2;
        if i == 0 || i == M
            H_tilda = H_tilda + 0.5*dt*H_i;
        else
            H_tilda = H_tilda + dt*H_i;
        end
    end

	f_ = x'*H_tilda*x;	
end

function [g_, h_] = my_constraints(x, p0, p1, X0, X1)

	% This is a linear constraint, so it need not be in this function; 
	% just for illustration 
	g_ = 0;
	
	% This is a nonlinear equality constraint
    A = [1 0 0 0 0 0 0 0;
         0 0 0 0 1 0 0 0;
         1 1 1 1 0 0 0 0;
         0 0 0 0 1 1 1 1;
         0 tan(X0) 0 0 0 -1 0 0;
         0 tan(X1) 2*tan(X1) 3*tan(X1) 0 -1 -2 -3];
    b = [p0(1) p0(2) p1(1) p1(2) 0 0]';
	h_ = A*x - b;
	% Surface area includes sides and the areas of the two circular ends
	
end