% final stage
function compare_risk_final_stage()
    % === Model Setup ===
    par = model.setup();
    par = model.gen_grids(par);

    r = par.r;
    beta = par.beta;
    sigma = par.sigma;
    A = 1.10;  % Fixed expected return

    w_grid = linspace(0.1, 25, 300);
    nw = length(w_grid);
    risk_scenarios = { [1.0], [1.0, 0.7], [1.0, 0.7, 0.4] };

    figure;
    for s = 1:length(risk_scenarios)
        risk_levels = risk_scenarios{s};
        V = zeros(1, nw);
        VW = zeros(1, nw);
        VE = zeros(1, nw);
        policy_k = zeros(1, nw);
        risk_p = ones(1, nw);

        ufun = @(c) ((c > 0) .* (c.^(1 - sigma)) ./ (1 - sigma)) - 1e5 * (c <= 0);
        tol = 1e-4; diff = 1; iter = 0; max_iter = 500;

        while diff > tol && iter < max_iter
            V_old = V;
            Vfun = @(wq) interp1(w_grid, V_old, wq, 'linear', 'extrap');

            for i = 1:nw
                w = w_grid(i);
                income = 1;
                agrid = linspace(0, w + income, 20);
                cons = w + income - agrid;
                VW(i) = max(ufun(cons) + beta * Vfun((1 + r) * agrid));

                kgrid = linspace(0, w, 20);
                best_val = -inf;
                for k = kgrid
                    c = w - k;
                    if c <= 0, continue; end
                    for p = risk_levels
                        x = A - (1 - p) * A / p;
                        w_low = x * k;
                        w_high = (A - x * (1 - p)) * k / p;
                        EV = p * Vfun(w_high) + (1 - p) * Vfun(w_low);
                        valE = ufun(c) + beta * EV;
                        if valE > best_val
                            best_val = valE;
                            policy_k(i) = k;
                            risk_p(i) = p;
                        end
                    end
                end
                VE(i) = best_val;
                V(i) = max(VE(i), VW(i));
            end
            diff = max(abs(V - V_old));
            iter = iter + 1;
        end

        % === Plot ===
        subplot(length(risk_scenarios), 2, 2*(s-1)+1);
        plot(w_grid, policy_k, 'LineWidth', 2);
        xlabel('Wealth'); ylabel('Investment k(w)');
        title(sprintf('Investment (Risk Options: %s)', mat2str(risk_levels)));

        subplot(length(risk_scenarios), 2, 2*(s-1)+2);
        plot(w_grid, risk_p, 'LineWidth', 2);
        xlabel('Wealth'); ylabel('Survival Prob p(w)');
        title(sprintf('Risk Level Chosen (Options: %s)', mat2str(risk_levels)));
    end

    sgtitle('Final Stage Comparison: Entrepreneurial Risk Choices');
end
