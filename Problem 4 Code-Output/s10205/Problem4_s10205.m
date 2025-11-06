%% Problem 4 - Patient s10205 Analysis
fprintf('\n=== Problem 4: Patient s10205 Analysis ===\n');

% Clear workspace to avoid contamination
clear; clc;

% Define parameters
TMAX_HR = 12;
TMAX_SEC = TMAX_HR * 3600;
OUT_DIR = '/Users/carolli/desktop/ICM Physiome/Project1';

% Patient s10205 specific paths
patient_id = 's10205';
patient_name = 'Patient #10205';
PATH_ABP_actual = '/Users/carolli/desktop/ICM Physiome/Project1/s10205-2631-06-14-11-38_ABP.txt';
PATH_NUM_actual = '/Users/carolli/desktop/ICM Physiome/Project1/s10205-2631-06-14-11-38n.txt';

fprintf('Processing %s\n', patient_name);
fprintf('ABP file: %s\n', PATH_ABP_actual);

% Check files exist
if ~exist(PATH_ABP_actual, 'file') || ~exist(PATH_NUM_actual, 'file')
    error('Files not found for %s', patient_name);
end

% Define manual correlation function
calculate_corr = @(x, y) sum((x - mean(x)) .* (y - mean(y))) / sqrt(sum((x - mean(x)).^2) * sum((y - mean(y)).^2));

try
    % Process this patient's data
    fprintf('Loading ABP data for %s...\n', patient_name);
    raw_p2 = readmatrix(PATH_ABP_actual, 'FileType', 'text');
    t_p2 = raw_p2(:,1);
    abp_p2 = raw_p2(:,2);
    
    % Limit to first 12 hours
    mask12_p2 = (t_p2 <= TMAX_SEC);
    time12_p2 = t_p2(mask12_p2);
    ABP12_p2 = abp_p2(mask12_p2);
    
    % Process ABP data
    fprintf('Processing ABP features for %s...\n', patient_name);
    t_on_p2 = wabp(ABP12_p2);
    feat_p2 = abpfeature(ABP12_p2, t_on_p2);
    beatq_p2 = jSQI(feat_p2, t_on_p2(1:end-1), ABP12_p2);
    
    % Save temporary mat file
    time = time12_p2; ABP = ABP12_p2;
    t_on = t_on_p2; feat = feat_p2; beatq = beatq_p2;
    temp_mat_file = sprintf('temp_%s_12h.mat', patient_id);
    save(temp_mat_file, 'time', 'ABP', 't_on', 'feat', 'beatq');
    
    % Run algorithms
    fprintf('Running algorithms for %s...\n', patient_name);
    [co_uncal_lilj_p2, to_min_p2, ~, ~] = estimateCO_v2(temp_mat_file, 5, 1);
    [co_uncal_pp_p2, ~, ~, ~] = estimateCO_v2(temp_mat_file, 2, 1);
    [co_uncal_ci_p2, ~, ~, ~] = estimateCO_v2(temp_mat_file, 7, 1);
    
    % Load TD data and calibrate
    fprintf('Loading thermodilution data for %s...\n', patient_name);
    tbl_p2 = readtable(PATH_NUM_actual, 'FileType', 'text', 'HeaderLines', 2);
    time_num_p2 = tbl_p2{:,1};
    CO_TD_p2 = tbl_p2{:,end};
    
    % Find nonzero CO values
    nz_idx_p2 = find(CO_TD_p2 ~= 0 & ~isnan(CO_TD_p2));
    
    if ~isempty(nz_idx_p2)
        td_time_sec_p2 = time_num_p2(nz_idx_p2);
        td_CO_Lmin_p2 = CO_TD_p2(nz_idx_p2);
        
        % Restrict to first 12h
        mask12_td_p2 = td_time_sec_p2 <= TMAX_SEC;
        td_time_sec_12h_p2 = td_time_sec_p2(mask12_td_p2);
        td_CO_Lmin_12h_p2 = td_CO_Lmin_p2(mask12_td_p2);
        td_time_hr_12h_p2 = td_time_sec_12h_p2 / 3600;
        
        fprintf('Found %d thermodilution measurements for %s\n', length(td_CO_Lmin_12h_p2), patient_name);
        
        if ~isempty(td_CO_Lmin_12h_p2)
            % Calibration using first measurement (C2 method)
            TD0_p2 = td_CO_Lmin_p2(1);
            t0_sec_p2 = td_time_sec_p2(1);
            [~, i_near_p2] = min(abs(to_min_p2 - (t0_sec_p2/60)));
            
            k_cal_lilj_p2 = TD0_p2 / co_uncal_lilj_p2(i_near_p2);
            k_cal_pp_p2 = TD0_p2 / co_uncal_pp_p2(i_near_p2);
            k_cal_ci_p2 = TD0_p2 / co_uncal_ci_p2(i_near_p2);
            
            % Apply calibration
            co_cal_lilj_p2 = k_cal_lilj_p2 * co_uncal_lilj_p2;
            co_cal_pp_p2 = k_cal_pp_p2 * co_uncal_pp_p2;
            co_cal_ci_p2 = k_cal_ci_p2 * co_uncal_ci_p2;
            
            % Calculate performance metrics
            valid_td_mask = (td_time_sec_12h_p2/60 >= min(to_min_p2)) & (td_time_sec_12h_p2/60 <= max(to_min_p2));
            
            if sum(valid_td_mask) == 0
                fprintf('WARNING: No TD measurements overlap with ABP time range for %s\n', patient_name);
                return;
            end
            
            % Use only valid TD measurements
            td_time_valid = td_time_sec_12h_p2(valid_td_mask);
            td_CO_valid = td_CO_Lmin_12h_p2(valid_td_mask);
            
            % Interpolate algorithm estimates at TD measurement times
            co_lilj_at_td_p2 = interp1(to_min_p2, co_cal_lilj_p2, td_time_valid/60, 'linear', 'extrap');
            co_pp_at_td_p2 = interp1(to_min_p2, co_cal_pp_p2, td_time_valid/60, 'linear', 'extrap');
            co_ci_at_td_p2 = interp1(to_min_p2, co_cal_ci_p2, td_time_valid/60, 'linear', 'extrap');
            
            % Calculate errors and performance metrics
            if ~any(isnan(co_lilj_at_td_p2)) && ~any(isnan(co_pp_at_td_p2)) && ~any(isnan(co_ci_at_td_p2))
                error_lilj_p2 = co_lilj_at_td_p2 - td_CO_valid;
                error_pp_p2 = co_pp_at_td_p2 - td_CO_valid;
                error_ci_p2 = co_ci_at_td_p2 - td_CO_valid;
                
                % RMS Error
                rms_lilj_p2 = sqrt(mean(error_lilj_p2.^2));
                rms_pp_p2 = sqrt(mean(error_pp_p2.^2));
                rms_ci_p2 = sqrt(mean(error_ci_p2.^2));
                
                % Correlation coefficients (manual calculation)
                corr_lilj_p2 = calculate_corr(co_lilj_at_td_p2, td_CO_valid);
                corr_pp_p2 = calculate_corr(co_pp_at_td_p2, td_CO_valid);
                corr_ci_p2 = calculate_corr(co_ci_at_td_p2, td_CO_valid);
                
                % Bland-Altman analysis (95% limits of agreement)
                % Liljestrand
                mean_diff_lilj = mean(error_lilj_p2);
                std_diff_lilj = std(error_lilj_p2);
                loa_lower_lilj = mean_diff_lilj - 1.96 * std_diff_lilj;
                loa_upper_lilj = mean_diff_lilj + 1.96 * std_diff_lilj;
                
                % Pulse Pressure
                mean_diff_pp = mean(error_pp_p2);
                std_diff_pp = std(error_pp_p2);
                loa_lower_pp = mean_diff_pp - 1.96 * std_diff_pp;
                loa_upper_pp = mean_diff_pp + 1.96 * std_diff_pp;
                
                % Corrected Impedance
                mean_diff_ci = mean(error_ci_p2);
                std_diff_ci = std(error_ci_p2);
                loa_lower_ci = mean_diff_ci - 1.96 * std_diff_ci;
                loa_upper_ci = mean_diff_ci + 1.96 * std_diff_ci;
                
                % Directional agreement (only if we have enough measurements)
                if length(td_CO_valid) >= 3
                    % Calculate directional changes in TD measurements
                    td_changes = diff(td_CO_valid);
                    td_directions = sign(td_changes);
                    
                    % Calculate corresponding changes in algorithm estimates
                    lilj_changes = diff(co_lilj_at_td_p2);
                    pp_changes = diff(co_pp_at_td_p2);
                    ci_changes = diff(co_ci_at_td_p2);
                    
                    lilj_directions = sign(lilj_changes);
                    pp_directions = sign(pp_changes);
                    ci_directions = sign(ci_changes);
                    
                    % Calculate directional agreement (exclude zero changes)
                    nonzero_mask = td_directions ~= 0;
                    if sum(nonzero_mask) > 0
                        dir_agree_lilj = 100 * sum(lilj_directions(nonzero_mask) == td_directions(nonzero_mask)) / sum(nonzero_mask);
                        dir_agree_pp = 100 * sum(pp_directions(nonzero_mask) == td_directions(nonzero_mask)) / sum(nonzero_mask);
                        dir_agree_ci = 100 * sum(ci_directions(nonzero_mask) == td_directions(nonzero_mask)) / sum(nonzero_mask);
                        n_directional = sum(nonzero_mask);
                    else
                        dir_agree_lilj = NaN;
                        dir_agree_pp = NaN;
                        dir_agree_ci = NaN;
                        n_directional = 0;
                    end
                else
                    dir_agree_lilj = NaN;
                    dir_agree_pp = NaN;
                    dir_agree_ci = NaN;
                    n_directional = 0;
                end
                
            else
                fprintf('WARNING: NaN values detected in interpolation for %s\n', patient_name);
                rms_lilj_p2 = NaN; rms_pp_p2 = NaN; rms_ci_p2 = NaN;
                corr_lilj_p2 = NaN; corr_pp_p2 = NaN; corr_ci_p2 = NaN;
                loa_lower_lilj = NaN; loa_upper_lilj = NaN;
                loa_lower_pp = NaN; loa_upper_pp = NaN;
                loa_lower_ci = NaN; loa_upper_ci = NaN;
                dir_agree_lilj = NaN; dir_agree_pp = NaN; dir_agree_ci = NaN;
                n_directional = 0;
            end
            
            % Display detailed performance comparison table
            fprintf('\n');
            fprintf('=================================================================\n');
            fprintf('Table: Agreement between thermodilution CO and CO-from-ABP algorithms\n');
            fprintf('Patient #10205 - First 12 hours (%d measurements)\n', length(td_CO_valid));
            fprintf('=================================================================\n');
            fprintf('\n');
            fprintf('                              95%% Limits of    RMS      Correlation  Directional\n');
            fprintf('Algorithm                     Agreement        Error    Coefficient  Agreement\n');
            fprintf('                              (±L/min)         (L/min)              (%%)\n');
            fprintf('-------------------------------------------------------------------------\n');
            
            if ~isnan(rms_lilj_p2)
                if ~isnan(dir_agree_lilj)
                    fprintf('Liljestrand                   %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
                            loa_lower_lilj, loa_upper_lilj, rms_lilj_p2, corr_lilj_p2, dir_agree_lilj);
                    fprintf('Pulse Pressure                %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
                            loa_lower_pp, loa_upper_pp, rms_pp_p2, corr_pp_p2, dir_agree_pp);
                    fprintf('Corrected Impedance           %+.2f/%+.2f      %.3f    %.3f        %.0f\n', ...
                            loa_lower_ci, loa_upper_ci, rms_ci_p2, corr_ci_p2, dir_agree_ci);
                else
                    fprintf('Liljestrand                   %+.2f/%+.2f      %.3f    %.3f        N/A\n', ...
                            loa_lower_lilj, loa_upper_lilj, rms_lilj_p2, corr_lilj_p2);
                    fprintf('Pulse Pressure                %+.2f/%+.2f      %.3f    %.3f        N/A\n', ...
                            loa_lower_pp, loa_upper_pp, rms_pp_p2, corr_pp_p2);
                    fprintf('Corrected Impedance           %+.2f/%+.2f      %.3f    %.3f        N/A\n', ...
                            loa_lower_ci, loa_upper_ci, rms_ci_p2, corr_ci_p2);
                end
            else
                fprintf('Unable to calculate metrics due to data issues\n');
            end
            
            fprintf('-------------------------------------------------------------------------\n');
            fprintf('\n');
            
            % Summary statistics
            fprintf('Summary Statistics:\n');
            fprintf('- Mean TD CO: %.2f ± %.2f L/min\n', mean(td_CO_valid), std(td_CO_valid));
            fprintf('- TD CO Range: %.2f - %.2f L/min\n', min(td_CO_valid), max(td_CO_valid));
            fprintf('- Number of directional changes analyzed: %d\n', n_directional);
            fprintf('\n');
            
            % Determine best algorithm
            if ~isnan(rms_lilj_p2)
                [~, best_rms_idx] = min([rms_lilj_p2, rms_pp_p2, rms_ci_p2]);
                [~, best_corr_idx] = max([corr_lilj_p2, corr_pp_p2, corr_ci_p2]);
                
                alg_names = {'Liljestrand', 'Pulse Pressure', 'Corrected Impedance'};
                fprintf('Best Performance:\n');
                fprintf('- Lowest RMS Error: %s (%.3f L/min)\n', alg_names{best_rms_idx}, min([rms_lilj_p2, rms_pp_p2, rms_ci_p2]));
                fprintf('- Highest Correlation: %s (%.3f)\n', alg_names{best_corr_idx}, max([corr_lilj_p2, corr_pp_p2, corr_ci_p2]));
                
                if ~isnan(dir_agree_lilj)
                    [~, best_dir_idx] = max([dir_agree_lilj, dir_agree_pp, dir_agree_ci]);
                    fprintf('- Best Directional Agreement: %s (%.0f%%)\n', alg_names{best_dir_idx}, max([dir_agree_lilj, dir_agree_pp, dir_agree_ci]));
                end
                fprintf('\n');
            end
            
            % Create plots
            fprintf('Creating plots for %s...\n', patient_name);
            to_hr_p2 = to_min_p2 / 60;
            
            % Plot 1: Liljestrand Algorithm
            figure('Color','w','Name',sprintf('%s Liljestrand', patient_name));
            hold on;
            plot(to_hr_p2, co_cal_lilj_p2, 'b-', 'LineWidth', 2);
            if ~isempty(td_time_hr_12h_p2)
                stem(td_time_hr_12h_p2, td_CO_Lmin_12h_p2, 'k', 'filled', 'LineWidth', 2);
            end
            ylabel('CO [L/min]');
            xlabel('Time [hours]');
            xlim([0 TMAX_HR]);
            grid on;
            title(sprintf('Liljestrand Algorithm (#5) - %s', patient_name));
            legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
            saveas(gcf, fullfile(OUT_DIR, sprintf('Problem4_%s_Liljestrand.png', patient_id)));
            
            % Plot 2: Pulse Pressure Algorithm
            figure('Color','w','Name',sprintf('%s Pulse Pressure', patient_name));
            hold on;
            plot(to_hr_p2, co_cal_pp_p2, 'r-', 'LineWidth', 2);
            if ~isempty(td_time_hr_12h_p2)
                stem(td_time_hr_12h_p2, td_CO_Lmin_12h_p2, 'k', 'filled', 'LineWidth', 2);
            end
            ylabel('CO [L/min]');
            xlabel('Time [hours]');
            xlim([0 TMAX_HR]);
            grid on;
            title(sprintf('Pulse Pressure Algorithm (#2) - %s', patient_name));
            legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
            saveas(gcf, fullfile(OUT_DIR, sprintf('Problem4_%s_PulsePressure.png', patient_id)));
            
            % Plot 3: Corrected Impedance Algorithm
            figure('Color','w','Name',sprintf('%s Corrected Impedance', patient_name));
            hold on;
            plot(to_hr_p2, co_cal_ci_p2, 'g-', 'LineWidth', 2);
            if ~isempty(td_time_hr_12h_p2)
                stem(td_time_hr_12h_p2, td_CO_Lmin_12h_p2, 'k', 'filled', 'LineWidth', 2);
            end
            ylabel('CO [L/min]');
            xlabel('Time [hours]');
            xlim([0 TMAX_HR]);
            grid on;
            title(sprintf('Corrected Impedance Algorithm (#7) - %s', patient_name));
            legend('Estimated CO', 'Thermodilution CO', 'Location', 'best');
            saveas(gcf, fullfile(OUT_DIR, sprintf('Problem4_%s_CorrectedImpedance.png', patient_id)));
            
            % Plot 4: Combined Comparison
            figure('Color','w','Name',sprintf('%s Comparison', patient_name));
            hold on;
            plot(to_hr_p2, co_cal_lilj_p2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Liljestrand (#5)');
            plot(to_hr_p2, co_cal_pp_p2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Pulse Pressure (#2)');
            plot(to_hr_p2, co_cal_ci_p2, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Corrected Impedance (#7)');
            if ~isempty(td_time_hr_12h_p2)
                stem(td_time_hr_12h_p2, td_CO_Lmin_12h_p2, 'k', 'filled', 'LineWidth', 2, ...
                     'DisplayName', 'Thermodilution CO');
            end
            
            ylabel('CO [L/min]');
            xlabel('Time [hours]');
            xlim([0 TMAX_HR]);
            grid on;
            legend('Location', 'best');
            title(sprintf('Algorithm Comparison - %s (First 12h)', patient_name));
            saveas(gcf, fullfile(OUT_DIR, sprintf('Problem4_%s_comparison.png', patient_id)));
            fprintf('Saved all figures for %s\n', patient_name);
            
        else
            fprintf('No thermodilution data in first 12h for %s\n', patient_name);
        end
    else
        fprintf('No thermodilution data found for %s\n', patient_name);
    end
    
    % Clean up temporary file
    delete(temp_mat_file);
    
catch ME
    fprintf('Error processing %s: %s\n', patient_name, ME.message);
    if exist('temp_mat_file', 'var') && exist(temp_mat_file, 'file')
        delete(temp_mat_file);
    end
end

fprintf('\n=== Patient %s Analysis Completed ===\n', patient_name);