% final stage
% Parameters
sigma = 2; beta = 0.958; phi = 1;
A = 1.10; r = 0.04;
w_min = 0.1; w_max = 25;
nw = 300; % fewer grid points
w_grid = linspace(w_min, w_max, nw);

% Utility
u = @(c) (c.^(1 - sigma)) / (1 - sigma);
u = @(c) (c > 0) .* u(c) - 1e6 * (c <= 0);  % penalize negative consumption

% Initialize
V = zeros(1, nw);
VE = zeros(1, nw);
VW = zeros(1, nw);
policy_k = zeros(1, nw);
risk_p = ones(1, nw);  % default to safe

% Iteration setup
tol = 1e-4; max_iter = 500;
diff = 1; iter = 0;

while diff > tol && iter < max_iter
    V_old = V;
    Vfun = @(wq) interp1(w_grid, V_old, wq, 'linear', 'extrap');

    for i = 1:nw
        w = w_grid(i);

        % ---- Worker value ----
        aw = linspace(0, w + phi, 20);
        cw = w + phi - aw;
        Vw1 = u(cw) + beta * Vfun((1 + r) * aw);
        [VW(i), ~] = max(Vw1);

        % ---- Entrepreneur value ----
        kw = linspace(0, w, 20);
        best_val = -inf;
        for k = kw
            c = w - k;
            if c <= 0, continue; end

            % Try 3 risk levels only
            for p = [1, 0.7, 0.4]
                x = A - (1 - p) * (A / p); % get x for given p
                wl = x * k;
                wh = (A - x * (1 - p)) * k / p;

                EV = p * Vfun(wh) + (1 - p) * Vfun(wl);
                val = u(c) + beta * EV;

                if val > best_val
                    best_val = val;
                    policy_k(i) = k;
                    risk_p(i) = p;
                end
            end
        end
        VE(i) = best_val;

        % Combine
        V(i) = max(VE(i), VW(i));
    end

    diff = max(abs(V - V_old));
    iter = iter + 1;
end

% ---- Plot results ----
figure;
subplot(2,1,1);
plot(w_grid, policy_k, 'LineWidth', 2);
xlabel('Wealth'); ylabel('Investment k(w)');
title('Entrepreneurial Investment Policy');

subplot(2,1,2);
plot(w_grid, risk_p, 'LineWidth', 2);
xlabel('Wealth'); ylabel('Survival Probability p(w)');
title('Risk Level Chosen by Entrepreneurs');
