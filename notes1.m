% first stage
% solve_stochastic_worker.m
clear; clc;

%% PARAMETERS
beta = 0.95;           % Discount factor
r = 0.03;              % Interest rate
sigma = 2.0;           % CRRA risk aversion

% Asset grid
grid_a = linspace(0.01, 20, 200)';
na = length(grid_a);

% Income states and transition matrix
y_vals = [1.0, 2.0];          % Low and high income
ny = length(y_vals);
Pi = [0.9 0.1;                % Transition: y1 -> [y1, y2]
      0.2 0.8];               % Transition: y2 -> [y1, y2]

% Utility function
utility = @(c) (c.^(1 - sigma)) / (1 - sigma);
if sigma == 1
    utility = @(c) log(c);
end

%% INITIALIZATION
V = zeros(na, ny);            % Value function
policy_c = zeros(na, ny);     % Consumption policy
policy_a_prime = zeros(na, ny);  % Savings policy

max_iter = 1000;
tol = 1e-4;

%% VALUE FUNCTION ITERATION
for iter = 1:max_iter
    V_new = zeros(na, ny);
    
    for iy = 1:ny
        y = y_vals(iy);
        
        for ia = 1:na
            a = grid_a(ia);
            a_next = grid_a;  % candidate next-period asset levels
            
            % Calculate consumption and utility
            c = (1 + r) * a + y - a_next;
            u_val = utility(c);
            u_val(c <= 0) = -1e10;  % penalize infeasible consumption
            
            % Expected continuation value
            EV = Pi(iy, :) * V';  % 1 x ny * ny x na = 1 x na
            total = u_val + beta * EV';
            
            [V_new(ia, iy), idx] = max(total);
            policy_c(ia, iy) = c(idx);
            policy_a_prime(ia, iy) = a_next(idx);
        end
    end
    
    % Convergence check
    if max(abs(V(:) - V_new(:))) < tol
        break;
    end
    V = V_new;
end

disp(['Converged in ', num2str(iter), ' iterations.']);

%% PLOTS: Policy Functions

figure;
for iy = 1:ny
    subplot(1, ny, iy);
    plot(grid_a, policy_c(:, iy), 'b-', 'LineWidth', 2); hold on;
    plot(grid_a, policy_a_prime(:, iy), 'r--', 'LineWidth', 2);
    title(['y = ', num2str(y_vals(iy))]);
    xlabel('Assets');
    legend('Consumption', 'Next-period Assets');
    grid on;
end
sgtitle('Stochastic Consumption-Savings Policy Functions');
