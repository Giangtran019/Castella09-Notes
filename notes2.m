%% second stage
%% Parameters
beta = 0.96;       % Discount factor
r = 0.04;          % Risk-free rate
A = 1.10;          % Expected project return
w_min = 0.01;      
w_max = 20;        
nw = 200;          % Grid points for wealth
sigma = 1e-5;      % Tolerance for convergence
max_iter = 1000;

w_grid = linspace(w_min, w_max, nw);   % Wealth grid

% Utility function (log utility)
u = @(c) log(c);
u_neg_inf = -1e10;

% Stochastic returns: risky project with two-point support
p = 0.5;             % Probability of low return
Z_low = 0.5;         % Low return
Z_high = 1.7;        % High return

%% Initialize
V = zeros(nw,1);           % Value function
policy_c = zeros(nw,1);    % Consumption policy
policy_k = zeros(nw,1);    % Capital investment policy

%% Value function iteration
for iter = 1:max_iter
    V_new = zeros(nw,1);
    for i = 1:nw
        w = w_grid(i);
        % Try all feasible investment levels
        obj = @(k) - (u(w - k) + beta * ...
            (p * interp1(w_grid, V, Z_low * k, 'linear', 0) + ...
             (1 - p) * interp1(w_grid, V, Z_high * k, 'linear', 0)));
        k_opt = fminbnd(obj, 0, w);
        c_opt = w - k_opt;
        
        V_new(i) = -obj(k_opt);
        policy_c(i) = c_opt;
        policy_k(i) = k_opt;
    end
    
    if max(abs(V_new - V)) < sigma
        disp(['Converged in ', num2str(iter), ' iterations']);
        break;
    end
    V = V_new;
end

%% Plot policy functions
figure;
plot(w_grid, policy_c, 'b-', 'LineWidth', 2); hold on;
plot(w_grid, policy_k, 'r--', 'LineWidth', 2);
xlabel('Wealth'); ylabel('Policy');
legend('Consumption', 'Investment in Project');
title('Policy Functions for Entrepreneur with Stochastic Return');
grid on;
