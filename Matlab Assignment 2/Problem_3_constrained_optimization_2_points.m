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

p0 = [0 0];
p1 = [10 -5];
X0 = 10*pi/180;
X1 = 50*pi/180;



%% Initial Guess

A_init = [1 p0(1) p0(1)^2 p0(1)^3;
          1 p1(1) p1(1)^2 p1(1)^3;
          0 1 2*p0(1) 3*p0(1)^2;
          0 1 2*p1(1) 3*p1(1)^2];
b_init = [p0(2) p1(2) tan(X0) tan(X1)]';

a = pinv(A_init)*b_init;
ax0 = p0(1);
ax1 = p1(1) - p0(1);
ay0 = a(1) + a(2)*ax0 + a(3)*ax0^2 + a(4)*ax0^3;
ay1 = a(2)*ax1 + 2*a(3)*ax0*ax1 + 3*a(4)*ax0^2*ax1;
ay2 = a(3)*ax1^2 + 3*a(4)*ax0*ax1^2;
ay3 = a(4)*ax1^3;

x_init = [ax0; ax1; 0; 0; ay0; ay1; ay2; ay3];
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
[x_opt, fval, exit_flag, solver_output] = ....
	fmincon(@(x) my_cost(x, 20, []), ...
	x_init, [], [], [], [], [0; 0], [], ...									% <--- The two zeros here are lower bounds for the decision variables
	@(x) my_constraints(x, p0, p1, X0, X1), ...
	solver_options);

fprintf('\n----- Optimal solution ----- \n')
fprintf('\t Optimal path coeff\t = %6.3f \n ', x_opt)
%fprintf('\t Height \t\t = %6.3f cm \n ', x_opt(2))


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