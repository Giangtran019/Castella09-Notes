% extension
% simulate_insurance_transitions.m
clear; clc;

%% Load saved policy functions and values
load('insurance_policy_data.mat');

%% GLOBAL PARAMETERS (add if not in .mat file)
r = 0.04;               % interest rate (for worker saving)
T = 100;                % number of periods
N = 1000;               % number of agents
nw = length(w_grid);

%% Initialize simulation arrays
wealth = zeros(T, N);
consumption = zeros(T, N);
occupation = zeros(T, N);  % 0 = worker, 1 = entrepreneur (uninsured), 2 = entrepreneur (insured)

% Initial wealth
wealth(1, :) = 1;

% Helper function to find nearest wealth index
find_index = @(x) max(1, sum(w_grid <= x));

%% Simulation loop
for t = 1:T-1
    for i = 1:N
        w = wealth(t, i);
        iw = find_index(w);

        % Compare occupational values
        val_worker = V_worker(iw);
        val_E = V_E(iw);
        val_EI = V_EI(iw);

        % Choose occupation
        [~, choice] = max([val_worker, val_E, val_EI]);
        occupation(t, i) = choice - 1;

        if choice == 1  % Worker
            c = policy_c_worker(iw);
            w_next = (1 + r) * (w - c);

        elseif choice == 2  % Entrepreneur (Uninsured)
            c = policy_c_E(iw);
            k = policy_k_E(iw);
            z = Z_L * (rand < p) + Z_H * (rand >= p);  % stochastic return
            w_next = z * k;

        else  % Entrepreneur (Insured)
            c = policy_c_EI(iw);
            k = policy_k_EI(iw);
            z = Z_insured * (rand < p) + Z_H * (rand >= p);
            w_next = z * k;
        end

        % Store and update
        consumption(t, i) = c;
        wealth(t+1, i) = max(w_next, 0.01);  % enforce minimum wealth
    end
end

% Final period consumption
for i = 1:N
    w = wealth(T, i);
    iw = find_index(w);
    val_worker = V_worker(iw);
    val_E = V_E(iw);
    val_EI = V_EI(iw);
    [~, choice] = max([val_worker, val_E, val_EI]);

    if choice == 1
        consumption(T, i) = policy_c_worker(iw);
    elseif choice == 2
        consumption(T, i) = policy_c_E(iw);
    else
        consumption(T, i) = policy_c_EI(iw);
    end
end

%% Plot results
avg_wealth = mean(wealth, 2);
avg_consumption = mean(consumption, 2);
share_worker = mean(occupation == 0, 2);
share_entrep = mean(occupation == 1, 2);
share_insured = mean(occupation == 2, 2);

figure;
subplot(3,1,1);
plot(1:T, avg_wealth, 'b', 'LineWidth', 2);
ylabel('Avg. Wealth'); title('Dynamic Transitions with Insurance'); grid on;

subplot(3,1,2);
plot(1:T, avg_consumption, 'r', 'LineWidth', 2);
ylabel('Avg. Consumption'); grid on;

subplot(3,1,3);
plot(1:T, share_worker, 'k', 'LineWidth', 2); hold on;
plot(1:T, share_entrep, 'r--', 'LineWidth', 2);
plot(1:T, share_insured, 'b:', 'LineWidth', 2);
legend('Worker', 'Entrepreneur', 'Entrepreneur (Insured)');
xlabel('Time'); ylabel('Occupational Share'); ylim([0 1]); grid on;
