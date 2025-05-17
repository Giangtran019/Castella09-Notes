%% second stage
% Simulation parameters
T = 100;                  % Number of periods
N = 1000;                 % Number of agents
wealth_path = zeros(T, N);
cons_path = zeros(T, N);
invest_path = zeros(T, N);

% Initial wealth
wealth_path(1,:) = 1;

% Simulate forward
for t = 1:T-1
    for n = 1:N
        w = wealth_path(t, n);
        if w <= w_min
            w = w_min;
        elseif w >= w_max
            w = w_max;
        end

        % Nearest wealth index
        w_idx = find(w_grid <= w, 1, 'last');
        if isempty(w_idx)
            w_idx = 1;
        end

        % Get policy decisions
        c = policy_c(w_idx);
        k = policy_k(w_idx);

        cons_path(t, n) = c;
        invest_path(t, n) = k;

        % Draw project return
        z = Z_low * (rand < p) + Z_high * (rand >= p);

        % Update next period wealth
        wealth_path(t+1, n) = z * k;
    end
end

%% Plot average wealth and consumption over time
mean_wealth = mean(wealth_path, 2);
mean_consumption = mean(cons_path, 2);

figure;
plot(1:T, mean_wealth, 'b-', 'LineWidth', 2); hold on;
plot(1:T, mean_consumption, 'r--', 'LineWidth', 2);
xlabel('Time');
ylabel('Mean Value');
legend('Wealth', 'Consumption');
title('Simulated Mean Wealth and Consumption Paths (Entrepreneur)');
grid on;
