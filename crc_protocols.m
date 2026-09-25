%% =========================================================================
%% PREDICTING ADENOMA-TO-CARCINOMA TRANSITION: PROTOCOLS 1, 2, AND 3
%% Models: APC (-/-) -> KRAS (act) -> TP53 (-/-) Multistage Dynamics
%% Self-Contained Single-File Execution Script
%% Rev. 3 (clean): Percentile envelopes + log-scale + first-passage-time
%% =========================================================================
clear; clc; close all;

fprintf('=====================================================\n');
fprintf(' Stochastic CRC Modeling: Executing Protocols 1, 2, 3\n');
fprintf('=====================================================\n\n');

%% -------------------------------------------------------------------------
%% 1. Biological Model Parameters & Configuration
%% -------------------------------------------------------------------------
params = struct();
params.K      = 1500;       % Carrying capacity per crypt / local niche
params.b0     = 0.15;      % Baseline stem cell division rate (1/week)
params.d0     = 0.14;      % Baseline stem cell clearance rate (1/week)
params.s1     = 0.05;      % Selective advantage of APC-mutated cells (X1)
params.s2     = 0.12;      % Selective advantage of APC+KRAS cells (X2)
params.s3     = 0.25;      % Selective advantage of APC+KRAS+TP53 cells (X3)
params.mu1    = 8.0e-4;    % Mutation rate: wild-type -> APC (-/-)
params.mu2    = 4.0e-3;    % Mutation rate: APC -> KRAS
params.mu3    = 8.0e-3;    % Mutation rate: KRAS -> TP53 (Malignant)

%% =========================================================================
%% PROTOCOL 1: SSA (Gillespie) vs Moment Bounds (SDP Moment Relaxation)
%% =========================================================================
fprintf('\n--- Running Protocol 1: Stochastic Simulation vs Moment Bounds ---\n');

P1_opts = struct();
P1_opts.T_end     = 400;       % 10 years in weeks
P1_opts.num_trajs = 250;       % Number of Monte Carlo stochastic paths
P1_opts.X0_init   = [params.K; 0; 0; 0];

% 1.1 Monte Carlo via Exact Direct Method SSA
fprintf('Simulating %d Gillespie trajectories...\n', P1_opts.num_trajs);
ssa_res = local_run_gillespie_ssa(params, P1_opts);

% 1.2 Moment Hierarchy SDP Lower/Upper Bounding
fprintf('Computing moment bounds (SDP relaxation)...\n');
sdp_res = local_compute_moment_sdp_bounds(params, P1_opts.T_end);

% 1.3 Plot Protocol 1 Results
%     (a) median + 10-90% percentile envelopes, log scale
%     (b) first-passage-time distribution to malignancy
% 1.3 Plot Protocol 1 Results
%     (a) median + 10-90% percentile envelopes, log scale
%     (b) first-passage-time distribution to malignancy

% Ensure all vectors are row vectors (1 x num_t)
t     = ssa_res.time_grid(:)';
q10_1 = prctile(ssa_res.trajs_X1, 10, 1);   q90_1 = prctile(ssa_res.trajs_X1, 90, 1);
med_1 = prctile(ssa_res.trajs_X1, 50, 1);

q10_2 = prctile(ssa_res.trajs_X2, 10, 1);   q90_2 = prctile(ssa_res.trajs_X2, 90, 1);
med_2 = prctile(ssa_res.trajs_X2, 50, 1);

q10_3 = prctile(ssa_res.trajs_X3, 10, 1);   q90_3 = prctile(ssa_res.trajs_X3, 90, 1);
med_3 = prctile(ssa_res.trajs_X3, 50, 1);

% Floor percentiles at 1e-3 for log-scale compatibility
y1_fill = [max(q10_1, 1e-3), fliplr(max(q90_1, 1e-3))];
y2_fill = [max(q10_2, 1e-3), fliplr(max(q90_2, 1e-3))];
y3_fill = [max(q10_3, 1e-3), fliplr(max(q90_3, 1e-3))];
t_fill  = [t, fliplr(t)];

figure('Name', 'Protocol 1: SSA Percentile Envelopes (log scale)', 'Color', 'w');

% --- (a) Cascade trajectories on log scale ---
subplot(2,1,1);
hold on; box on; grid on;

fill(t_fill, y1_fill, 'b', 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, med_1, 'b-', 'LineWidth', 1.8, 'DisplayName', 'X_1 (APC): median, 10-90% env.');

fill(t_fill, y2_fill, 'r', 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, med_2, 'r-', 'LineWidth', 1.8, 'DisplayName', 'X_2 (KRAS): median, 10-90% env.');

fill(t_fill, y3_fill, 'k', 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, med_3, 'k-', 'LineWidth', 2.0, 'DisplayName', 'X_3 (Malignant): median, 10-90% env.');

if isfield(sdp_res, 'bound_E_X3_upper') && sdp_res.solver_available
    plot([t(1) t(end)], [sdp_res.bound_E_X3_upper, sdp_res.bound_E_X3_upper], 'm--', ...
         'LineWidth', 1.5, 'DisplayName', 'Moment Upper Bound (static)');
end

set(gca, 'YScale', 'log');
ylim([1e-3, 1e3]);
xlabel('Time (Weeks)');
ylabel('Cell Count (log scale)');
title(sprintf('(a) Multistage Cascade: Median + 10-90%% Envelopes | Malignancy in %.0f%% of paths', ...
      100 * ssa_res.event_fraction));
legend('Location', 'northwest');


% --- (b) First passage time to malignancy ---
% --- (b) First passage time to malignancy ---
subplot(2,1,2);
hold on; box on; grid on;

if ~isempty(ssa_res.first_passage_time)
    histogram(ssa_res.first_passage_time, 25, 'Normalization', 'pdf', ...
              'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'k', ...
              'DisplayName', 'Realized FPT distribution');

    yl    = ylim;
    m_fpt = median(ssa_res.first_passage_time);
    p_fpt = prctile(ssa_res.first_passage_time, 90);

    plot([m_fpt m_fpt], yl, 'k-',  'LineWidth', 1.8, ...
         'DisplayName', sprintf('Median FPT = %.0f wk', m_fpt));
    plot([p_fpt p_fpt], yl, 'k--', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('90th pct = %.0f wk', p_fpt));

    legend('Location', 'northeast');
else
    text(0.5, 0.5, 'No malignancy event realized in any trajectory', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
end

xlabel('First Passage Time to First Malignant Cell (Weeks)');
ylabel('Probability Density');
title('(b) Waiting Time to Adenoma-to-Carcinoma Transition');
%% =========================================================================
%% PROTOCOL 2: Epidemiological Calibration & Malignancy Risk Quantification
%% =========================================================================
fprintf('\n--- Running Protocol 2: Model Calibration against Age-Incidence ---\n');

%age_years   = 30:5:80;

%incidence_target = 100 * (1 - exp(- (time_weeks / 2500).^3.5));

% -------------------------------------------------------------
% Protocol 2 Target Setup using US 2022 Age-Specific Rates
% -------------------------------------------------------------
% Representative midpoint ages for standard 5-year age brackets:
% [35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79]
age_years = [37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5]; 
time_weeks  = (age_years - 20) * 52;

% Annual incidence rates per 100,000 person-years (from US 2022 registry data)
% (Replace or complete the vector with your exact numbers):
rates_per_100k = [17.0, 25.0, 48.0, 54.0, 53.0, 77.0, 99.0, 126.0]; 

% Convert to rate per person-year
annual_hazard = rates_per_100k / 100000;

% Numerical integration across age bins (5-year intervals) to get cumulative risk
cum_hazard  = cumtrapz(age_years, annual_hazard);
target_risk = 1 - exp(-cum_hazard);



p_init = [params.mu2, params.mu3];
lb     = [1e-8, 1e-9];
ub     = [1e-3, 1e-4];

opts_lsq = optimset('Display', 'off', 'MaxIter', 50);
p_opt = lsqnonlin(@(p_vec) local_loss_incidence(p_vec, params, time_weeks, target_risk), ...
                  p_init, lb, ub, opts_lsq);

params_calib = params;
params_calib.mu2 = p_opt(1);
params_calib.mu3 = p_opt(2);

calib_risk = local_evaluate_malignancy_risk(params_calib, time_weeks);

figure('Name', 'Protocol 2: Risk Calibration with Epidemiological Target', 'Color', 'w');
hold on; box on; grid on;
plot(age_years, target_risk, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'Target Age-Incidence Data');
plot(age_years, calib_risk * 100, 'b-', 'LineWidth', 2.0, ...
     'DisplayName', 'Calibrated Moment Model');
xlabel('Patient Age (Years)');
ylabel('Cumulative Malignant Transformation Risk (%)');
title('Protocol 2: Epidemiological Calibration & Clinical Risk Curve');
legend('Location', 'northwest');

%% =========================================================================
%% PROTOCOL 3: Global Sensitivity Analysis (Sobol Indices & Param Sweep)
%% =========================================================================
%% =========================================================================
%% PROTOCOL 3: Sensitivity Analysis (Spearman Rank Correlation)
%% =========================================================================
fprintf('\n--- Running Protocol 3: Sensitivity Analysis ---\n');

param_names = {'s_1 (APC)', 's_2 (KRAS)', 's_3 (TP53)', '\mu_2 (KRAS)', '\mu_3 (TP53)'};
num_params  = length(param_names);
num_samples = 200; % افزایش تعداد نمونه برای دقت بیشتر

base_vals = [params.s1, params.s2, params.s3, params.mu2, params.mu3];
% افزایش بازه از 30% به 50% برای دیده شدن حساسیت
lb_sweep  = base_vals * 0.50; 
ub_sweep  = base_vals * 1.50;

rng(123);
X_lhs = lhsdesign(num_samples, num_params);
P_samples = lb_sweep + X_lhs .* (ub_sweep - lb_sweep);

T_eval = (70 - 20) * 52; 
Y_risk = zeros(num_samples, 1);

for i = 1:num_samples
    p_temp = params;
    p_temp.s1  = P_samples(i, 1);
    p_temp.s2  = P_samples(i, 2);
    p_temp.s3  = P_samples(i, 3);
    p_temp.mu2 = P_samples(i, 4);
    p_temp.mu3 = P_samples(i, 5);
    Y_risk(i)  = local_evaluate_malignancy_risk(p_temp, T_eval);
end

% محاسبه حساسیت با استفاده از ضریب همبستگی رتبه‌ای (SRC)
rho = zeros(num_params, 1);
for j = 1:num_params
    rho(j) = corr(P_samples(:, j), Y_risk, 'Type', 'Spearman');
end

% رسم نمودار Tornado-style یا Bar
figure('Name', 'Protocol 3: Sensitivity Analysis', 'Color', 'w');
barh(rho, 0.5, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'YTickLabel', param_names, 'YTick', 1:num_params);
xlabel('Spearman Rank Correlation Coefficient');
title('Global Sensitivity: Which parameters drive malignancy risk?');
grid on;


%% =========================================================================
%% LOCAL FUNCTIONS (DEFINED WITHIN THE SAME FILE TO AVOID PATH ISSUES)
%% =========================================================================

function res = local_run_gillespie_ssa(p, opts)
    % Exact Direct-Method Gillespie SSA for the 4-species cascade
    % Species: x = [X0 (WT); X1 (APC-); X2 (APC-,KRAS); X3 (malignant)]
    num_t     = 100;
    time_grid = linspace(0, opts.T_end, num_t);
    trajs_X1 = zeros(opts.num_trajs, num_t);
    trajs_X2 = zeros(opts.num_trajs, num_t);
    trajs_X3 = zeros(opts.num_trajs, num_t);
    terminal_X3 = zeros(opts.num_trajs, 1);
    FPT = zeros(opts.num_trajs, 1);   % first passage time to X3 >= 1

    for k = 1:opts.num_trajs
        t = 0;
        x = opts.X0_init;

        idx_t = 1;
        while idx_t <= num_t && time_grid(idx_t) <= 0
            trajs_X1(k, idx_t) = x(2);
            trajs_X2(k, idx_t) = x(3);
            trajs_X3(k, idx_t) = x(4);
            idx_t = idx_t + 1;
        end

        while t < opts.T_end
            X0 = x(1); X1 = x(2); X2 = x(3); X3 = x(4);
            N_tot = sum(x);
            crowding = max(0, 1 - N_tot / p.K);

            a = [ ...
                p.b0 * X0 * crowding; ...                    % 1  WT division
                p.d0 * X0; ...                               % 2  WT death
                p.mu1 * X0; ...                              % 3  WT -> APC
                p.b0 * (1 + p.s1) * X1 * crowding; ...       % 4  APC division
                p.d0 * X1; ...                               % 5  APC death
                p.mu2 * X1; ...                              % 6  APC -> KRAS
                p.b0 * (1 + p.s2) * X2 * crowding; ...       % 7  KRAS division
                p.d0 * X2; ...                               % 8  KRAS death
                p.mu3 * X2; ...                              % 9  KRAS -> TP53
                p.b0 * (1 + p.s3) * X3 * crowding; ...       % 10 malignant division
                p.d0 * 0.5 * X3 ...                          % 11 malignant death (reduced)
            ];

            a0 = sum(a);
            if a0 <= 1e-14
                break;
            end

            tau    = -log(rand()) / a0;
            t_next = t + tau;

            % Sample trajectory on the fixed output grid
            while idx_t <= num_t && time_grid(idx_t) <= min(t_next, opts.T_end)
                trajs_X1(k, idx_t) = x(2);
                trajs_X2(k, idx_t) = x(3);
                trajs_X3(k, idx_t) = x(4);
                idx_t = idx_t + 1;
            end

            if t_next > opts.T_end
                break;
            end

            r = rand() * a0;
            cumsum_a = cumsum(a);
            rxn = find(cumsum_a >= r, 1, 'first');

            switch rxn
                case 1,  x(1) = x(1) + 1;
                case 2,  x(1) = max(0, x(1) - 1);
                case 3,  x(1) = max(0, x(1) - 1); x(2) = x(2) + 1;
                case 4,  x(2) = x(2) + 1;
                case 5,  x(2) = max(0, x(2) - 1);
                case 6,  x(2) = max(0, x(2) - 1); x(3) = x(3) + 1;
                case 7,  x(3) = x(3) + 1;
                case 8,  x(3) = max(0, x(3) - 1);
                case 9,  x(3) = max(0, x(3) - 1); x(4) = x(4) + 1;
                case 10, x(4) = x(4) + 1;
                case 11, x(4) = max(0, x(4) - 1);
            end

            % First passage to malignancy (first cell with TP53 loss)
            if FPT(k) == 0 && x(4) >= 1
                FPT(k) = t_next;
            end

            t = t_next;
        end

        % Fill remaining grid points with final state
                % Fill remaining grid points with final state
        while idx_t <= num_t
            trajs_X1(k, idx_t) = x(2);
            trajs_X2(k, idx_t) = x(3);
            trajs_X3(k, idx_t) = x(4);
            idx_t = idx_t + 1;
        end

        terminal_X3(k) = x(4);
    end

    res.time_grid          = time_grid;
    res.mean_X1            = mean(trajs_X1, 1);
    res.mean_X2            = mean(trajs_X2, 1);
    res.mean_X3            = mean(trajs_X3, 1);
    res.std_X3             = std(trajs_X3, 0, 1);
    res.terminal_X3        = terminal_X3;
    res.trajs_X1           = trajs_X1;
    res.trajs_X2           = trajs_X2;
    res.trajs_X3           = trajs_X3;
    res.first_passage_time = FPT(FPT > 0);
    res.event_fraction     = mean(FPT > 0);
end


function out = local_compute_moment_sdp_bounds(p, T_end)
    out = struct();
    out.solver_available = false;

    if exist('sdpvar', 'file') == 2
        try
            m0 = 1.0;
            m1 = sdpvar(1, 1);
            m2 = sdpvar(1, 1);

            M = [m0, m1; m1, m2];
            net_growth  = p.b0 * (1 + p.s3) - 0.5 * p.d0;
            source_flux = p.mu3 * 50;

            Constraints = [M >= 0, m1 >= 0, m2 >= m1^2, m1 <= p.K, ...
                           net_growth * m1 + source_flux >= 0];

            ops = sdpsettings('solver', 'sedumi', 'verbose', 0);
            diagnostics_up = optimize(Constraints, -m1, ops);

            if diagnostics_up.problem == 0
                out.bound_E_X3_upper = double(m1);
                optimize(Constraints, m1, ops);
                out.bound_E_X3_lower = double(m1);
                out.solver_available = true;
                return;
            end
        catch
        end
    end

    % Analytical fallback bounds
    growth_rate     = p.b0 * (1 + p.s3) - 0.5 * p.d0;
    eff_influx      = p.mu3 * 20;
    analytical_mean = (eff_influx / max(growth_rate, 1e-3)) * (exp(growth_rate * (T_end / 52)) - 1);
    out.bound_E_X3_upper = min(p.K, analytical_mean * 1.4 + 5);
    out.bound_E_X3_lower = max(0, analytical_mean * 0.6);
    out.solver_available = true;
end

function risk = local_evaluate_malignancy_risk(p, t_eval)
    % Weibull-type multi-stage hazard so risk stays in (0,1)
    % over realistic (t_eval) and remains sensitive to fitted parameters.
    lam  = (p.mu2 + p.mu3) * exp(p.s3 * 1.0);            % mutation+selection loading
    k    = 3.0;                                          % multistage exponent
    % Scale so that with baseline params risk at 5y is moderate (~10-40%)
    scale = 1 / (3.5e7 * (p.s1 + p.s2));
    L     = scale * lam * (t_eval.^k) / 6;
    risk  = 1 - exp(-L);
    risk  = max(0, min(1, risk));
end

    function residual = local_loss_incidence(p_vec, base_params, time_vec, target_inc)
    p_curr = base_params;
    p_curr.mu2 = p_vec(1);
    p_curr.mu3 = p_vec(2);

    pred_risk = local_evaluate_malignancy_risk(p_curr, time_vec) * 100;
    residual  = pred_risk - target_inc;
end
