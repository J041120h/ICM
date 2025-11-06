
%% -------------------------- PROBLEM 5 ---------------------------------
%% Continuous CO estimation using Parlikar (2007) estimator (#14)
%
% This section can be APPENDED to your existing general code (do not modify earlier parts).
% It repeats the minimal preprocessing to save a 0–12h segment and then calls
% estimateCO_v2 with estID = 14. You may set PATH_ABP / PATH_NUM as desired.
%
% Inputs to provide before running this section:
%   PATH_ABP : path to the patient's ABP waveform file (two columns: t [s], ABP [mmHg])
%   PATH_NUM : path to the patient's numeric file (thermodilution column last)
%   OUT_DIR  : output directory for figures/results

fprintf('\n=== Starting Problem 5: Continuous CO estimation (Parlikar #14) ===\n');

if ~exist('PATH_ABP','var') || ~exist('PATH_NUM','var') || ~exist('OUT_DIR','var')
    error(['Please define PATH_ABP, PATH_NUM, and OUT_DIR before running Problem 5.\n' ...
           'Example:\n  PATH_ABP = ''/path/to/s00020-..._ABP.txt'';\n' ...
           '  PATH_NUM = ''/path/to/s00020-...n.txt'';\n' ...
           '  OUT_DIR  = ''/desired/output/dir'';']);
end
if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

% -------------------- DEFINE BASIC PARAMETERS --------------------------
Fs = 125;
TMAX_HR = 12;
TMAX_SEC = TMAX_HR * 3600;

% ------------------------------- READ DATA -----------------------------
raw = readmatrix(PATH_ABP,'FileType','text'); % [time(s), ABP(mmHg)]
t   = raw(:,1);
abp = raw(:,2);
t   = t(:);
abp = abp(:);

% ------------------------- TRIM ABP (0–12h) ----------------------------
mask12 = (t <= TMAX_SEC);
time12 = t(mask12);
ABP12  = abp(mask12);

% ---------------------- RUN 2analyze FUNCTIONS -------------------------
fprintf('Running WABP / ABPFEATURE / JSQL for 12-hour data (Q5)\n');
t_on  = wabp(ABP12);                      % beat onsets
feat  = abpfeature(ABP12,t_on);           % beat features
beatq = jSQI(feat,t_on(1:end-1),ABP12);   % SQI

% Save features for estimator v2
time = time12; ABP = ABP12;
mat_out = fullfile(OUT_DIR,'patient_first12h_Q5.mat');
save(mat_out,'time','ABP','t_on','feat','beatq');

% ---------------------- ESTIMATE CO (UNCALIBRATED) ---------------------
fprintf('Running Parlikar estimator (#14)\n');
[co_uncal,to_min,~,fea] = estimateCO_v2(mat_out,14,1);
to_hr = to_min / 60;

% Extract features for plotting
Psys   = fea(:,2);
Pdias  = fea(:,4);
PP     = fea(:,5);
MAP    = fea(:,6);
Period = fea(:,7);
HR     = 60 * Fs ./ Period;

% ---------------------- LOAD THERMODILUTION CO -------------------------
fprintf('Reading numeric file for CO measurements (Q5)\n');
tbl = readtable(PATH_NUM,'FileType','text','HeaderLines',2);
time_num_all = tbl{:,1};              % seconds
CO_TD_all    = tbl{:,end};            % last column = CO [L/min]

% find nonzero CO values
nz_idx_all     = find(CO_TD_all ~= 0 & ~isnan(CO_TD_all));
td_time_sec_all = time_num_all(nz_idx_all);
td_CO_Lmin_all  = CO_TD_all(nz_idx_all);

% restrict to first 12h
mask12_TD       = td_time_sec_all <= TMAX_SEC;
td_time_sec_12h = td_time_sec_all(mask12_TD);
td_CO_Lmin_12h  = td_CO_Lmin_all(mask12_TD);
td_time_hr_12h  = td_time_sec_12h / 3600;

% ---------------------- CALIBRATION C2 / C3 ----------------------------
if isempty(td_time_sec_all)
    warning('No nonzero thermodilution CO found in numeric file. Using uncalibrated scale.');
    k_cal = 1;
else
    TD0    = td_CO_Lmin_all(1);                 % first TD value
    t0_sec = td_time_sec_all(1);                % time of first TD
    [~,i_near] = min(abs(to_min - (t0_sec/60)));
    est0 = co_uncal(i_near);
    k_cal = TD0 / est0;                         % calibration factor
    fprintf('Calibration factor k = %.3f (first TD at %.2f hr)\n',k_cal,t0_sec/3600);
end
co_cal = k_cal * co_uncal;

% ------------------------------ PLOTS ----------------------------------
f = figure('Color','w');
subplot(2,1,1);
plot(to_hr, co_cal, 'LineWidth', 1.1); hold on;
scatter(td_time_hr_12h, td_CO_Lmin_12h, 30, 'filled'); grid on; box on;
xlabel('Time (hr)'); ylabel('CO (L/min)'); 
legend('Parlikar (#14) calibrated','Thermodilution','Location','best');
title('Continuous CO (Parlikar 2007) vs Thermodilution');

subplot(2,1,2);
plot(to_hr, PP, 'LineWidth', 1); hold on;
plot(to_hr, MAP, 'LineWidth', 1);
plot(to_hr, HR,  'LineWidth', 1);
grid on; box on;
xlabel('Time (hr)'); ylabel('Feature value');
legend('PP','MAP','HR','Location','best');
title('ABP Features (first 12 h)');

print(fullfile(OUT_DIR,'Problem5_Parlikar_CO_and_Features.png'),'-dpng','-r300');
close(gcf);

fprintf('Saved Q5 results to %s\n', OUT_DIR);
