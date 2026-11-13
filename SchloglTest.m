clear; clc; close all;

%% 1. Parameters (Strictly Bistable Regime with Real Bimodality)
params.k1 = 200;      % 0 -> X
params.k2 = 3.5;      % X -> 0
params.k3 = 0.03;     % 2X -> 3X
params.k4 = 0.0001;   % 3X -> 2X

opts.T_end     = 80;     % زمان کافی برای واگرایی مسیرها
opts.dt        = 0.2;
opts.num_trajs = 500;
opts.X0_init   = 250;    % شروع دقیقاً بر روی پشته ناپایدار (Unstable Separatrix)

%% 2. Run Gillespie SSA
fprintf('Running Gillespie SSA for Schlögl system (Bistable)...\n');
res = run_schlogl_gillespie(params, opts);

%% 3. Plot Dynamics & Bistability (Publication Ready)
figure('Color', 'w', 'Position', [100, 100, 750, 480]);
t = res.time_grid;

% 10th-90th Percentile Band (Shows Bimodality / Branching)
fill([t, fliplr(t)], [res.percentiles(1,:), fliplr(res.percentiles(3,:))], ...
     [0.85 0.90 0.95], 'EdgeColor', 'none', 'DisplayName', '10th–90th Percentile');
hold on;

% Median and Mean
plot(t, res.percentiles(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Median SSA');
plot(t, res.mean_X, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean SSA');

xlabel('Time ($t$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Molecule Count ($X$)', 'Interpreter', 'latex', 'FontSize', 12);
title('Schlögl Model: Bimodal Stochastic Switching', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 10);
grid on; box on;
xlim([0 opts.T_end]);
ylim([0 700]);

%% =========================================================================
%% Local Helper Function (Gillespie SSA)
%% =========================================================================
function res = run_schlogl_gillespie(p, opts)
    t_grid = 0:opts.dt:opts.T_end;
    N_steps = length(t_grid);
    trajs = zeros(opts.num_trajs, N_steps);
    
    for i = 1:opts.num_trajs
        t = 0;
        x = opts.X0_init;
        idx = 1;
        trajs(i, idx) = x;
        
        while t < opts.T_end && idx < N_steps
            % Propensities for Schlögl system
            a1 = p.k1;
            a2 = p.k2 * x;
            a3 = p.k3 * max(0, x * (x - 1) / 2);
            a4 = p.k4 * max(0, x * (x - 1) * (x - 2) / 6);
            
            a0 = a1 + a2 + a3 + a4;
            if a0 <= 0, break; end
            
            tau = -log(rand()) / a0;
            t = t + tau;
            
            while idx < N_steps && t_grid(idx + 1) <= t
                idx = idx + 1;
                trajs(i, idx) = x;
            end
            
            if idx >= N_steps, break; end
            
            % Reaction selection
            r = rand() * a0;
            if r <= a1
                x = x + 1;                  % 0 -> X
            elseif r <= a1 + a2
                x = x - 1;                  % X -> 0
            elseif r <= a1 + a2 + a3
                x = x + 1;                  % 2X -> 3X
            else
                x = x - 1;                  % 3X -> 2X
            end
        end
        trajs(i, idx:end) = x;
    end
    
    res.time_grid = t_grid;
    res.mean_X = mean(trajs, 1);
    res.percentiles = prctile(trajs, [10, 50, 90], 1);
end
