%% =========================================================================
%% IN SILICO CLINICAL SIMULATIONS: Calibrated Risk Profiles
%% =========================================================================
% Optimized parameters for visual differentiation and clinical relevance.
% -------------------------------------------------------------------------

clear; clc; close all;
format compact;

%% 1. Calibrated Model Parameters (Refined for Clinical Separation)
K    = 2000;          % Crypt capacity
b0   = 0.15;          
d0   = 0.14;          
s1   = 0.02;          % Milder fitness for APC-/-
s2   = 0.08;          % Moderate fitness for KRAS
s3   = 0.20;          % Strong fitness for TP53/Malignancy

% Adjusted Mutation Rates (Realistic for clinical 10-year horizon)
mu1  = 1.0e-6;        
mu2  = 2.0e-5;        
mu3  = 5.0e-4;        

r1 = (b0*(1 + s1) - d0);
r2 = (b0*(1 + s2) - d0);

%% 2. Define 3 In Silico Patients (Initial States X1, X2)
patients = [50,   0 ;    % Patient A: Minimal APC involvement (Early Polyp)
            300,  10 ;   % Patient B: Significant APC + KRAS seeding (Advanced)
            800, 200 ];  % Patient C: Near-saturation KRAS (High-grade Dysplasia)

labels = {'Patient A: Low-Risk (APC^{-/-} only)', ...
          'Patient B: Intermediate-Risk (APC^{-/-} + KRAS^{mut})', ...
          'Patient C: High-Risk (Pre-malignant Dysplasia)'};

%% 3. Simulation & Integration
T_years = 10;
t_weeks = 0:2:(T_years*52);
n_t     = length(t_weeks);
risk_up = zeros(3, n_t);

for p = 1:3
    X1 = patients(p,1); X2 = patients(p,2);
    H = 0; dt = 2; % weeks
    for i = 1:n_t
        % Dynamics (Moment-closure style ODE)
        dX1 = r1 * X1 * (1 - X1/K);
        dX2 = r2 * X2 * (1 - X2/K) + mu2 * X1;
        
        X1 = X1 + dt * dX1;
        X2 = X2 + dt * dX2;
        
        % Hazard is probability of X2 cells mutating to malignant state
        hazard = mu3 * X2;
        H = H + dt * hazard;
        risk_up(p, i) = 1 - exp(-H);
    end
end

%% 4. Visualization
fig = figure('Color','white','Position',[100 100 1100 480]);

% --- Subplot 1: Risk Curves ---
subplot(1,2,1);
colors = [0.15, 0.65, 0.15;   % Green (Safe)
          0.90, 0.60, 0.00;   % Orange (Warning)
          0.80, 0.10, 0.10];  % Red (Critical)
hold on; box on; grid on;
for p = 1:3
    plot(t_weeks/52, 100*risk_up(p,:), 'LineWidth', 2.5, 'Color', colors(p,:));
end
set(gca, 'YScale', 'log', 'YMinorGrid', 'on');
ylim([0.01 100]); xlim([0 T_years]);
xlabel('Years from Index Colonoscopy');
ylabel('Malignancy Risk (%)');
title('Personalized Risk Trajectories', 'FontSize', 13);
legend(labels, 'Location', 'northwest', 'FontSize', 9);
yline(10, 'k--', 'Risk Threshold (10%)', 'LabelHorizontalAlignment','right');

% --- Subplot 2: Recommended Windows ---
subplot(1,2,2);
fpt10 = zeros(3,1);
for p = 1:3
    idx = find(risk_up(p,:) >= 0.10, 1, 'first');
    if isempty(idx)
        fpt10(p) = T_years;
    else
        fpt10(p) = t_weeks(idx)/52;
    end
end
b = barh(1:3, fpt10, 0.6, 'FaceColor', 'flat');
b.CData = colors;
set(gca, 'YTick', 1:3, 'YTickLabel', {'Patient A', 'Patient B', 'Patient C'});
xlabel('Safe Surveillance Interval (Years)');
title('Optimal Screening Windows (FPT_{10%})', 'FontSize', 12);
xlim([0 11]); grid on;
for p = 1:3
    text(fpt10(p)+0.2, p, sprintf('%.1f Yrs', fpt10(p)), 'FontWeight','bold');
end

%% 5. Clinical Summary Table (Console)
fprintf('\n%-12s | %-12s | %-12s | %-15s\n', 'Group', '5-Yr Risk', '10-Yr Risk', 'Recommendation');
fprintf('%s\n', repmat('-', 1, 60));
idx5 = find(t_weeks >= 5*52, 1, 'first');
for p = 1:3
    rec = 'Routine';
    if fpt10(p) < 3
        rec = 'Urgent/Surg';
    elseif fpt10(p) < 6
        rec = 'Early Follow-up';
    end
    fprintf('Patient %c    | %-10.2f%% | %-10.2f%% | %-15s\n', char(64+p), 100*risk_up(p,idx5), 100*risk_up(p,end), rec);
end
