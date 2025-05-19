% stage 1
% simulate_stochastic_worker.m
clear;

%% LOAD POLICIES FROM SOLVER
load('policy_data.mat', 'policy_a_prime', 'policy_c', 'grid_a', 'y_vals', 'Pi');

%% SIMULATION PARAMETERS
T = 200;                  % Number of periods
N = 1000;                 % Number of agents

% Initialize arrays
assets = zeros(T, N);          % asset holdings
consumption = zeros(T, N);     % consumption
income_state = zeros(T, N);    % income state index (1 or 2)
income_level = zeros(T, N);    % actual income

% Initial conditions
a0 = 1;                          % start with 1 unit of assets
assets(1, :) = a0;
income_state(1, :) = 1;          % start with low income
income_level(1, :) = y_vals(1);

% Helper function to find the nearest index on the asset grid
find_index = @(x, grid) max(1, sum(grid <= x));

%% MAIN SIMULATION LOOP
for t = 1:T-1
    for n = 1:N
        a = assets(t, n);
        iy = income_state(t, n);

        % Locate nearest asset index
        ia = find_index(a, grid_a);

        % Record consumption
        consumption(t, n) = policy_c(ia, iy);

        % Income state transition
        p = rand;
        cdf = cumsum(Pi(iy, :));
        iy_next = find(p <= cdf, 1);
        income_state(t+1, n) = iy_next;
        income_level(t+1, n) = y_vals(iy_next);

        % Update asset using policy function
        assets(t+1, n) = policy_a_prime(ia, iy);
    end
end

% Final period consumption
for n = 1:N
    a = assets(T, n);
    iy = income_state(T, n);
    ia = find_index(a, grid_a);
    consumption(T, n) = policy_c(ia, iy);
end

%% AVERAGE PLOTS

mean_assets = mean(assets, 2);
mean_consumption = mean(consumption, 2);

figure;
subplot(2,1,1);
plot(1:T, mean_assets, 'b', 'LineWidth', 2);
xlabel('Time'); ylabel('Average Assets');
title('Mean Asset Path over 200 Periods');
grid on;

subplot(2,1,2);
plot(1:T, mean_consumption, 'r', 'LineWidth', 2);
xlabel('Time'); ylabel('Average Consumption');
title('Mean Consumption Path over 200 Periods');
grid on;