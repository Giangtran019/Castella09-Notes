% extension
% solve_insured_vs_uninsured.m
clear; clc;

%% Parameters
beta = 0.96;
r = 0.04;
sigma = 1.0;  % Log utility

% Wealth grid
w_min = 0.01; w_max = 20; nw = 200;
w_grid = linspace(w_min, w_max, nw);

% TUNED project return settings
p = 0.5;        % probability of low return
Z_L = 0.9;      % tuned upward from 0.5
Z_H = 1.2;      % tuned downward from 1.7

% Insurance policy
Z_insured = 1.0;    % tuned floor
premium = 0.02;

% Utility function
u = @(c) log(c);
penalty = -1e10;

%% Worker Value Function
V_worker = zeros(nw, 1);
policy_c_worker = zeros(nw, 1);

for iter = 1:500
    V_new = zeros(nw, 1);
    for i = 1:nw
        w = w_grid(i);
        c_try = linspace(0.001, w, nw);
        a_next = (1 + r) * (w - c_try);
        util = u(c_try); util(c_try <= 0) = penalty;
        V_interp = interp1(w_grid, V_worker, a_next, 'linear', 0);
        V_total = util + beta * V_interp;
        [V_new(i), idx] = max(V_total);
        policy_c_worker(i) = c_try(idx);
    end
    if max(abs(V_worker - V_new)) < 1e-4, break; end
    V_worker = V_new;
end

%% Entrepreneur – Uninsured
V_E = zeros(nw, 1);
policy_c_E = zeros(nw, 1);
policy_k_E = zeros(nw, 1);

for iter = 1:500
    V_new = zeros(nw, 1);
    for i = 1:nw
        w = w_grid(i);
        if w <= 0
            V_new(i) = penalty;
            continue;
        end
        obj = @(k) - (u(max(w - k, 1e-8)) + beta * ...
            (p * interp1(w_grid, V_E, max(Z_L * k, w_min), 'linear', 0) + ...
             (1 - p) * interp1(w_grid, V_E, max(Z_H * k, w_min), 'linear', 0)));
        k_opt = fminbnd(obj, 0, w);
        c_opt = max(w - k_opt, 1e-8);
        V_new(i) = -obj(k_opt);
        policy_c_E(i) = c_opt;
        policy_k_E(i) = k_opt;
    end
    if max(abs(V_E - V_new)) < 1e-4, break; end
    V_E = V_new;
end

%% Entrepreneur – Insured
V_EI = zeros(nw, 1);
policy_c_EI = zeros(nw, 1);
policy_k_EI = zeros(nw, 1);

for iter = 1:500
    V_new = zeros(nw, 1);
    for i = 1:nw
        w = w_grid(i);
        if w <= premium
            V_new(i) = penalty;
            continue;
        end
        obj = @(k) - (u(max(w - k - premium, 1e-8)) + beta * ...
            (p * interp1(w_grid, V_EI, max(Z_insured * k, w_min), 'linear', 0) + ...
             (1 - p) * interp1(w_grid, V_EI, max(Z_H * k, w_min), 'linear', 0)));
        k_upper = max(w - premium, 0);
        k_opt = fminbnd(obj, 0, k_upper);
        c_opt = max(w - k_opt - premium, 1e-8);
        V_new(i) = -obj(k_opt);
        policy_c_EI(i) = c_opt;
        policy_k_EI(i) = k_opt;
    end
    if max(abs(V_EI - V_new)) < 1e-4, break; end
    V_EI = V_new;
end

%% Save for simulation
save('insurance_policy_data.mat', ...
    'w_grid', 'V_worker', 'V_E', 'V_EI', ...
    'policy_c_worker', 'policy_c_E', 'policy_k_E', ...
    'policy_c_EI', 'policy_k_EI', ...
    'Z_L', 'Z_H', 'Z_insured', 'p', 'premium');

%% Plot: Occupational Choice With and Without Insurance
is_entrepreneur_uninsured = V_E > V_worker;
is_entrepreneur_insured = V_EI > V_worker;

figure;
hold on;
plot(w_grid, double(is_entrepreneur_uninsured), 'r--', 'LineWidth', 2, 'DisplayName', 'Uninsured');
plot(w_grid, double(is_entrepreneur_insured), 'b-', 'LineWidth', 2, 'DisplayName', 'Insured');
xlabel('Wealth'); ylabel('Entrepreneur Chosen (1=True)');
title('Occupational Choice With Tuned Risk and Insurance');
legend('show');
ylim([-0.1 1.1]); grid on;

fprintf('Share choosing uninsured entrepreneur: %.2f%%\n', 100 * mean(is_entrepreneur_uninsured));
fprintf('Share choosing insured entrepreneur: %.2f%%\n', 100 * mean(is_entrepreneur_insured));
