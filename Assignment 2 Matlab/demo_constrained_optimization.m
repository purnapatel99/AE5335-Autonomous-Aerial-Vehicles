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

surface_area	= 1000;														% cm^3, surface area should be exactly equal to this number
max_height		= 20;														% height should be no greater than this number

%% Initial Guess

x_init = [1;20];
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
	fmincon(@(x) my_cost(x, [], []), ...
	x_init, [], [], [], [], [0; 0], [], ...									% <--- The two zeros here are lower bounds for the decision variables
	@(x) my_constraints(x, surface_area, max_height), ...
	solver_options);

fprintf('\n----- Optimal solution ----- \n')
fprintf('\t Radius of base\t = %6.3f cm \n ', x_opt(1))
fprintf('\t Height \t\t = %6.3f cm \n ', x_opt(2))


function f_ = my_cost(x, parameter_1, parameter_2)
	% The variables 'parameter_1' etc are not used here, but included
	% merely to illustrate the syntax for passing parameters into the cost
	% function. By default the syntax for this cost function has only one
	% argument, namely, the decision variable x. There can be as many
	% parameters as needed passed as arguments to this function.
	
	f_ = -pi * x(1)^2 * x(2);	
end

function [g_, h_] = my_constraints(x, parameter_1, parameter_2)
	% This function returns nonlinear inequality constraints g_ and
	% nonlinear equality constraints h_

	% Let's say parameter_1 is the desired surface area, whereas
	% 'parameter_2' is the max height

	% This is a linear constraint, so it need not be in this function; 
	% just for illustration 
	g_ = x(2) - parameter_2;
	
	% This is a nonlinear equality constraint
	h_ = 2*pi*x(1)*x(2) + 2*pi*x(1)^2 - parameter_1;
	% Surface area includes sides and the areas of the two circular ends
	
end