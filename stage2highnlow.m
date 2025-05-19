% proposition1_thresholds.m
clear; clc;

%% PARAMETERS
r = 0.63;                     % interest rate (1 + r = 1.63)
A = 2.16;                     % project return
beta = 0.95;                  % discount factor
phi = 1.0;                    % labor income (ϕ)

% Safe utility function to avoid log(0)
u = @(c) (c > 0) .* log(c) + (c <= 0) * (-1e10);

% Wealth grid
w0_grid = linspace(0.01, 5, 300);

% Period-1 value function V1(w1)
V1 = @(w1) u(w1) + beta * u((1 + r) * w1);

% Initialize value functions
V0E = zeros(size(w0_grid));  % entrepreneur value at t = 0
V0W = zeros(size(w0_grid));  % worker value at t = 0

% Loop over wealth grid
for i = 1:length(w0_grid)
    w0 = w0_grid(i);

    %% ENTREPRENEUR VALUE FUNCTION
    % Decision variables: [k0, b0]
    obj_E = @(x) - ( u(w0 - x(1) - x(2)) + ...
                     beta * V1(A * x(1) + (1 + r) * x(2)) );

    % Constraints: k0 + b0 <= w0, k0 >= 0, b0 >= 0
    Aineq = [1 1]; bineq = w0;
    lb = [0 0];

    % Safe initial guess (to ensure consumption > 0)
    x0 = [w0 / 3, w0 / 3];

    % Run optimization
    [~, fval_E] = fmincon(obj_E, x0, Aineq, bineq, [], [], lb, [], [], ...
        optimoptions('fmincon', 'Display', 'off'));

    V0E(i) = -fval_E;

    %% WORKER VALUE FUNCTION
    obj_W = @(w1) - ( u(w0 + phi - w1 / (1 + r)) + beta * V1(w1) );
    [~, fval_W] = fminbnd(obj_W, 0.01, w0 + phi);
    V0W(i) = -fval_W;
end

%% FIND THRESHOLDS w0^L and w0^H
diff_V = V0E - V0W;

% Find lower threshold where entrepreneur value > worker
idx_L = find(diff_V > 0, 1, 'first');
idx_H = find(diff_V < 0 & (1:length(diff_V))' > idx_L, 1, 'first');

% Extract threshold values
w0_L = NaN; w0_H = NaN;
if ~isempty(idx_L), w0_L = w0_grid(idx_L); end
if ~isempty(idx_H), w0_H = w0_grid(idx_H); end

fprintf('Estimated Thresholds:\n');
fprintf('  w0_L = %.4f\n', w0_L);
fprintf('  w0_H = %.4f\n', w0_H);

%% PLOT VALUE FUNCTIONS
figure;
plot(w0_grid, V0E, 'r-', 'LineWidth', 2); hold on;
plot(w0_grid, V0W, 'b--', 'LineWidth', 2);

% Plot threshold markers
if ~isnan(w0_L), xline(w0_L, 'k--', 'LineWidth', 1.5); end
if ~isnan(w0_H), xline(w0_H, 'k--', 'LineWidth', 1.5); end

legend('Entrepreneur V_0^E', 'Worker V_0^W', 'w_0^L', 'w_0^H', ...
       'Location', 'best');
xlabel('Initial Wealth w_0');
ylabel('Value');
title('Proposition 1: Occupational Choice Thresholds at t = 0');
grid on;