%% -------------------------- PROBLEM 4 ----------------------------------
%% Compare Liljestrand with two other algorithms for Patient #20

fprintf('\n=== Starting Problem 4: Algorithm Comparison for Patient #20 ===\n');

% Define basic parameters
TMAX_HR = 12;
TMAX_SEC = TMAX_HR * 3600;

% Use the same 12-hour data from Problem 3
% The .mat file should already exist from Problem 3
if ~exist('patient20_first12h.mat', 'file')
    fprintf('Creating patient20_first12h.mat file...\n');
    
    % Recreate the 12-hour data processing from Problem 3
    PATH_ABP = '/Users/carolli/desktop/ICM Physiome/Project1/s00020-2567-03-30-17-47_ABP.txt';
    
    raw = readmatrix(PATH_ABP,'FileType','text');
    t = raw(:,1);
    abp = raw(:,2);
    
    % Trim to first 12 hours
    mask12 = (t <= TMAX_SEC);
    time = t(mask12);
    ABP = abp(mask12);
    
    % Process ABP data
    t_on = wabp(ABP);
    feat = abpfeature(ABP, t_on);
    beatq = jSQI(feat, t_on(1:end-1), ABP);
    
    % Save the .mat file
    save('patient20_first12h.mat', 'time', 'ABP', 't_on', 'feat', 'beatq');
    fprintf('Created patient20_first12h.mat successfully\n');
end

% ---------------------- LOAD THERMODILUTION DATA ------------------------
fprintf('Loading thermodilution data...\n');
PATH_NUM = '/Users/carolli/desktop/ICM Physiome/Project1/s00020-2567-03-30-17-47n.txt';
tbl = readtable(PATH_NUM,'FileType','text','HeaderLines',2);

time_num_all = tbl{:,1};  % seconds
CO_TD_all = tbl{:,end}; % last column = CO [L/min]

% find nonzero CO values
nz_idx_all = find(CO_TD_all ~= 0 & ~isnan(CO_TD_all));
td_time_sec_all = time_num_all(nz_idx_all); % timestamps in seconds
td_CO_Lmin_all  = CO_TD_all(nz_idx_all); % CO values [L/min]

% restrict to first 12h
mask12_TD = td_time_sec_all <= TMAX_SEC;
td_time_sec_12h = td_time_sec_all(mask12_TD);
td_CO_Lmin_12h  = td_CO_Lmin_all(mask12_TD);
td_time_hr_12h  = td_time_sec_12h / 3600;

fprintf('Found %d thermodilution measurements in first 12h\n', length(td_time_sec_12h));

% ---------------------- ESTIMATE CO WITH ALL THREE ALGORITHMS -----------
fprintf('Running all three algorithms...\n');

% Algorithm 5: Liljestrand
[co_uncal_lilj, to_min, ~, fea] = estimateCO_v2('patient20_first12h.mat', 5, 1);
to_hr = to_min / 60;

% Algorithm 2: Pulse Pressure (Windkessel)
[co_uncal_pp, ~, ~, ~] = estimateCO_v2('patient20_first12h.mat', 2, 1);

% Algorithm 7: Corrected Impedance (Wesseling)
[co_uncal_ci, ~, ~, ~] = estimateCO_v2('patient20_first12h.mat', 7, 1);

% ---------------------- APPLY INDIVIDUAL CALIBRATION (C2 METHOD) --------
if ~isempty(td_time_sec_all)
    TD0 = td_CO_Lmin_all(1); % first measured TCO value
    t0_sec = td_time_sec_all(1); % time of first measured CO
    [~, i_near] = min(abs(to_min - (t0_sec/60)));
    
    % Calculate SEPARATE calibration factors for each algorithm
    k_cal_lilj = TD0 / co_uncal_lilj(i_near);
    k_cal_pp = TD0 / co_uncal_pp(i_near);
    k_cal_ci = TD0 / co_uncal_ci(i_near);
    
    fprintf('Calibration factors: Lilj=%.3f, PP=%.6f, CI=%.9f\n', k_cal_lilj, k_cal_pp, k_cal_ci);
else
    k_cal_lilj = 1; k_cal_pp = 1; k_cal_ci = 1;
    fprintf('No thermodilution data found, using unit calibration\n');
end

% Apply INDIVIDUAL calibration to each algorithm
co_cal_lilj = k_cal_lilj * co_uncal_lilj;
co_cal_pp = k_cal_pp * co_uncal_pp;
co_cal_ci = k_cal_ci * co_uncal_ci;

% ---------------------- CALCULATE PERFORMANCE METRICS ------------------
if ~isempty(td_time_sec_12h)
    % Interpolate algorithm estimates at TD measurement times
    co_lilj_at_td = interp1(to_min, co_cal_lilj, td_time_sec_12h/60, 'linear');
    co_pp_at_td = interp1(to_min, co_cal_pp, td_time_sec_12h/60, 'linear');
    co_ci_at_td = interp1(to_min, co_cal_ci, td_time_sec_12h/60, 'linear');
    
    % Calculate errors
    error_lilj = co_lilj_at_td - td_CO_Lmin_12h;
    error_pp = co_pp_at_td - td_CO_Lmin_12h;
    error_ci = co_ci_at_td - td_CO_Lmin_12h;
    
    % Calculate performance metrics
    rms_lilj = sqrt(mean(error_lilj.^2));
    rms_pp = sqrt(mean(error_pp.^2));
    rms_ci = sqrt(mean(error_ci.^2));
    
    mae_lilj = mean(abs(error_lilj));
    mae_pp = mean(abs(error_pp));
    mae_ci = mean(abs(error_ci));
    
    % Calculate 95% limits of agreement (mean ± 1.96*std)
    bias_lilj = mean(error_lilj);
    bias_pp = mean(error_pp);
    bias_ci = mean(error_ci);
    
    loa_lower_lilj = bias_lilj - 1.96*std(error_lilj);
    loa_upper_lilj = bias_lilj + 1.96*std(error_lilj);
    
    loa_lower_pp = bias_pp - 1.96*std(error_pp);
    loa_upper_pp = bias_pp + 1.96*std(error_pp);
    
    loa_lower_ci = bias_ci - 1.96*std(error_ci);
    loa_upper_ci = bias_ci + 1.96*std(error_ci);
    
    % Manual correlation calculation
    x_lilj = co_lilj_at_td - mean(co_lilj_at_td);
    y_td = td_CO_Lmin_12h - mean(td_CO_Lmin_12h);
    corr_lilj = sum(x_lilj .* y_td) / sqrt(sum(x_lilj.^2) * sum(y_td.^2));
    
    x_pp = co_pp_at_td - mean(co_pp_at_td);
    corr_pp = sum(x_pp .* y_td) / sqrt(sum(x_pp.^2) * sum(y_td.^2));
    
    x_ci = co_ci_at_td - mean(co_ci_at_td);
    corr_ci = sum(x_ci .* y_td) / sqrt(sum(x_ci.^2) * sum(y_td.^2));
    
    % Calculate directional agreement (simplified version)
    % For each consecutive pair of TD measurements, check if algorithm follows same direction
    n_directional = 0;
    correct_lilj = 0; correct_pp = 0; correct_ci = 0;
    
    if length(td_CO_Lmin_12h) > 1
        for i = 1:length(td_CO_Lmin_12h)-1
            td_direction = sign(td_CO_Lmin_12h(i+1) - td_CO_Lmin_12h(i));
            lilj_direction = sign(co_lilj_at_td(i+1) - co_lilj_at_td(i));
            pp_direction = sign(co_pp_at_td(i+1) - co_pp_at_td(i));
            ci_direction = sign(co_ci_at_td(i+1) - co_ci_at_td(i));
            
            if td_direction ~= 0  % Only count non-zero changes
                n_directional = n_directional + 1;
                if lilj_direction == td_direction; correct_lilj = correct_lilj + 1; end
                if pp_direction == td_direction; correct_pp = correct_pp + 1; end
                if ci_direction == td_direction; correct_ci = correct_ci + 1; end
            end
        end
    end
    
    % Calculate directional agreement percentages
    if n_directional > 0
        dir_agree_lilj = (correct_lilj / n_directional) * 100;
        dir_agree_pp = (correct_pp / n_directional) * 100;
        dir_agree_ci = (correct_ci / n_directional) * 100;
    else
        dir_agree_lilj = 0; dir_agree_pp = 0; dir_agree_ci = 0;
    end
    
    % Display performance comparison table (similar to Table 3 format)
    fprintf('\n');
    fprintf('=================================================================\n');
    fprintf('Table: Agreement between thermodilution CO and CO-from-ABP algorithms\n');
    fprintf('Patient #20 - First 12 hours (%d measurements)\n', length(td_CO_Lmin_12h));
    fprintf('=================================================================\n');
    fprintf('\n');
    fprintf('                              95%% Limits of    RMS      Correlation  Directional\n');
    fprintf('Algorithm                     Agreement        Error    Coefficient  Agreement\n');
    fprintf('                              (±L/min)         (L/min)              (%%)\n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Liljestrand                   %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
            loa_lower_lilj, loa_upper_lilj, rms_lilj, corr_lilj, dir_agree_lilj);
    fprintf('Pulse Pressure                %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
            loa_lower_pp, loa_upper_pp, rms_pp, corr_pp, dir_agree_pp);
    fprintf('Corrected Impedance           %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
            loa_lower_ci, loa_upper_ci, rms_ci, corr_ci, dir_agree_ci);
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('\n');
    
    % Summary statistics
    fprintf('Summary Statistics:\n');
    fprintf('- Mean TD CO: %.2f ± %.2f L/min\n', mean(td_CO_Lmin_12h), std(td_CO_Lmin_12h));
    fprintf('- TD CO Range: %.2f - %.2f L/min\n', min(td_CO_Lmin_12h), max(td_CO_Lmin_12h));
    fprintf('- Number of directional changes analyzed: %d\n', n_directional);
    fprintf('\n');
    
    % Determine best algorithm
    [~, best_rms_idx] = min([rms_lilj, rms_pp, rms_ci]);
    [~, best_corr_idx] = max([corr_lilj, corr_pp, corr_ci]);
    [~, best_dir_idx] = max([dir_agree_lilj, dir_agree_pp, dir_agree_ci]);
    
    alg_names = {'Liljestrand', 'Pulse Pressure', 'Corrected Impedance'};
    fprintf('Best Performance:\n');
    fprintf('- Lowest RMS Error: %s (%.3f L/min)\n', alg_names{best_rms_idx}, min([rms_lilj, rms_pp, rms_ci]));
    fprintf('- Highest Correlation: %s (%.3f)\n', alg_names{best_corr_idx}, max([corr_lilj, corr_pp, corr_ci]));
    fprintf('- Best Directional Agreement: %s (%.0f%%)\n', alg_names{best_dir_idx}, max([dir_agree_lilj, dir_agree_pp, dir_agree_ci]));
    fprintf('\n');
end

% ---------------------- CREATE INDIVIDUAL ALGORITHM PLOTS --------------
fprintf('Creating individual algorithm plots...\n');

OUT_DIR = '/Users/carolli/desktop/ICM Physiome/Project1';

% Plot 1: Liljestrand Algorithm
figure('Color','w','Name','Liljestrand Algorithm');
hold on;
plot(to_hr, co_cal_lilj, 'b-', 'LineWidth', 2);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2);
end
ylabel('CO [L/min]');
xlabel('Time [hours]');
xlim([0 TMAX_HR]);
grid on;
title('Liljestrand Algorithm (#5) - Patient #20');
legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
saveas(gcf, fullfile(OUT_DIR, 'Problem4_Liljestrand.png'));

% Plot 2: Pulse Pressure Algorithm
figure('Color','w','Name','Pulse Pressure Algorithm');
hold on;
plot(to_hr, co_cal_pp, 'r-', 'LineWidth', 2);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2);
end
ylabel('CO [L/min]');
xlabel('Time [hours]');
xlim([0 TMAX_HR]);
grid on;
title('Pulse Pressure Algorithm (#2) - Patient #20');
legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
saveas(gcf, fullfile(OUT_DIR, 'Problem4_PulsePressure.png'));

% Plot 3: Corrected Impedance Algorithm
figure('Color','w','Name','Corrected Impedance Algorithm');
hold on;
plot(to_hr, co_cal_ci, 'g-', 'LineWidth', 2);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2);
end
ylabel('CO [L/min]');
xlabel('Time [hours]');
xlim([0 TMAX_HR]);
grid on;
title('Corrected Impedance Algorithm (#7) - Patient #20');
legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
saveas(gcf, fullfile(OUT_DIR, 'Problem4_CorrectedImpedance.png'));

% Plot 4: Combined Comparison
figure('Color','w','Name','Algorithm Comparison');
hold on;
plot(to_hr, co_cal_lilj, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Liljestrand (#5)');
plot(to_hr, co_cal_pp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Pulse Pressure (#2)');
plot(to_hr, co_cal_ci, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Corrected Impedance (#7)');
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2, ...
         'DisplayName', 'Thermodilution CO');
end
ylabel('CO [L/min]');
xlabel('Time [hours]');
xlim([0 TMAX_HR]);
grid on;
legend('Location', 'best');
title('Algorithm Comparison - Patient #20 (First 12h)');
saveas(gcf, fullfile(OUT_DIR, 'Problem4_Comparison.png'));

fprintf('Saved all figures to %s\n', OUT_DIR);

fprintf('\n=== Problem 4 Completed ===\n');