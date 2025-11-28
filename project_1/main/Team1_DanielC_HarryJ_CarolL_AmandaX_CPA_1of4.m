% Project 1 – Estimating Cardiac Output from ABP 
% Team 1: Daniel Chen, Harry Jiang, Carol Li, Amanda Xu
% Problem 1 & Problem 3 Combined Script


%% -------------------------- PROBLEM 1 ----------------------------------
% Reproduce traces similar to Figure 1 in paper for Patient #20
clear; clc;

fprintf('\n=== Starting Problem 1: Producing ABP traces for Patient #20 ===');

% --------------------------DEFINING VARIABLES------------------------
% define paths for inputs/outputs. EDIT before running
PATH_ABP = ['s00020-2567-03-30-17-47_ABP.txt'];   % ABP waveform file
PATH_NUM = 's00020-2567-03-30-17-47n.txt';      % numeric file for Problem 3
OUT_DIR  = '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/Q5'; % path for output

% create output folder if it does not exist
if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

% define time anchors from the ABP records
% as we plot the first 20 pulses star777ting at 10 h and 11 h, we find
% timepoint for 10hours and 11hours in the ABP file.
T10 = 10*3600;   % 10 hours
T11 = 11*3600;   % 11 hours

% -------------------------------READ DATA--------------------------------
% read the ABP file: col 1 is time in seconds, col 2 is ABP in mmHg
raw = readmatrix(PATH_ABP,'FileType','text');
t   = raw(:,1);
abp = raw(:,2);
t   = t(:); % enforce vectors to be column vectors
abp = abp(:);

% ---------------------------BEAT ONSET DETECTION-------------------------
% obtain onset time for each beat using helper function
OnsetTimes = wabp(abp);

% ------------------- SELECT 21 ONSETS (FOR 20 BEATS) ---------------------
% General Algo:
%   1) Find the first sample index in the time vector at or after the anchor.
%   2) Find the first detected onset at or after that sample index.
%   3) Take 21 consecutive onset indices to ensure exactly 20 beats.

% for 10 h
start10_idx = find(t >= T10, 1, 'first');
k10 = find(OnsetTimes >= start10_idx, 1, 'first');
onsets10 = OnsetTimes(k10:(k10+20));   % 21 onsets = 20 beats

% For 11 h
start11_idx = find(t >= T11, 1, 'first');
k11 = find(OnsetTimes >= start11_idx, 1, 'first');
onsets11 = OnsetTimes(k11:(k11+20));   % 21 onsets = 20 beats

% --------------- EXTRACTING FEATURES ---------------------------
% from abpfeature.m,
%   col  9: End of systole time by 0.3*sqrt(RR) method
%   col 11: End of systole time by lowest negative slope method
feat10 = abpfeature(abp, onsets10); 
feat11 = abpfeature(abp, onsets11);


% -------------------- CALCULATE SIGNAL QUALITY ------------------------
% calculate ABP waveform signal quality index using helper function
[BeatQ10, ~] = jSQI(feat10, onsets10(1:end-1), abp);
[BeatQ11, ~] = jSQI(feat11, onsets11(1:end-1), abp);

% --------------------------- PLOT (10 h) -------------------------------
fprintf('\nPlotting traces for 10h and 11h\n');

% reproduce the plot like figure 1
i0 = onsets10(1);     % first onset sample
i1 = onsets10(end);   % 21st onset sample
tt = t(i0:i1);
yy = abp(i0:i1);

on20   = onsets10(1:end-1);   % 20 onset markers (one per beat)
eos03  = feat10(:,9);   % end-systole by 0.3*sqrt(RR)
eosMin = feat10(:,11);  % end-systole by lowest negative slope method

f1 = figure('Color','w'); % create a new figure
plot(tt, yy, 'LineWidth', 1); hold on;  % ABP trace
plot(t(on20),abp(on20),'*','MarkerSize',5); % onsets
plot(t(eos03), abp(eos03), 'xr', 'MarkerSize', 7, 'LineWidth', 1.2); % 0.3*sqrt(RR) ends
plot(t(eosMin),abp(eosMin),'o','MarkerSize',5);  % lowest neg slope ends
xlabel('Time (s)');
ylabel('ABP (mmHg)');
title('Patient #20, Starting at 10 h');
box on; grid on;

% save as png
print(fullfile(OUT_DIR,'Problem1_10h.png'),'-dpng','-r300');
close(gcf);  

% --------------------------- PLOT (11 h) -------------------------------
% the logic is the exact same as above.
i0 = onsets11(1);
i1 = onsets11(end);
tt = t(i0:i1);
yy = abp(i0:i1);

on20   = onsets11(1:end-1);
eos03  = feat11(:,9);
eosMin = feat11(:,11);

f2 = figure('Color','w');
plot(tt, yy, 'LineWidth', 1); hold on;
plot(t(on20), abp(on20), '*', 'MarkerSize', 5);
%plot(t(eos03), abp(eos03), 'xr', 'MarkerSize', 7, 'LineWidth', 1.2);
plot(t(eos03), abp(eos03), 'xr', 'MarkerSize', 7, 'LineWidth', 1.2); % 0.3*sqrt(RR) ends
plot(t(eosMin), abp(eosMin), 'o', 'MarkerSize', 5);
xlabel('Time (s)');
ylabel('ABP (mmHg)');
title('Patient #20, Starting at 11 h');
box on; grid on;

ax = findall(gcf,'type','axes');
set(ax, 'Toolbar', []);
print(fullfile(OUT_DIR,'Problem1_11h.png'),'-dpng','-r300');
close(gcf);
fprintf('Saved Figure to %s\n',OUT_DIR);
fprintf('\n=== Problem 1 Completed ===\n');


%% -------------------------- PROBLEM 3 ---------------------------------
%% Continuous CO estimation using Liljestrand (#5)

fprintf('\n=== Starting Problem 3: Continuous CO estimation (Liljestrand) for Patient 20 ===\n');

% -------------------- DEFINE BASIC PARAMETERS --------------------------
Fs = 125;
TMAX_HR = 12;
TMAX_SEC = TMAX_HR * 3600;

% ------------------------- TRIM ABP (0–12h) -----------------------------
mask12 = (t <= TMAX_SEC); % Limit ABP signal to first 12 hours only
time12 = t(mask12);
ABP12  = abp(mask12);


fprintf('Running WABP / ABPFEATURE / JSQL for 12-hour data\n');
t_on  = wabp(ABP12); % Detect beat onsets in ABP waveform
feat  = abpfeature(ABP12,t_on); % Extract beat-by-beat ABP features
beatq = jSQI(feat,t_on(1:end-1),ABP12);  % Compute beat-level signal quality 

% save into a .mat file for use by the estimator
time = time12;
ABP  = ABP12;
save('patient20_first12h.mat','time','ABP','t_on','feat','beatq'); 


% ---------------------- ESTIMATE CO (UNCALIBRATED) ---------------------
fprintf('Running Liljestrand estimator (#5)\n');
[co_uncal,to_min,~,fea] = estimateCO_v2('patient20_first12h.mat',5,1);
to_hr = to_min / 60;

% Extract features for plotting
Psys   = fea(:,2);
Pdias  = fea(:,4);
PP     = fea(:,5);
MAP    = fea(:,6);
Period = fea(:,7);
HR     = 60 * Fs ./ Period;

% ---------------------- LOAD THERMODILUTION CO --------------------------
fprintf('Reading numeric file for CO measurements\n');
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

% ---------------------- CALIBRATION C2(Sun 2009)/C3 (Sun 2005) -------------------------------
if isempty(td_time_sec_all)
    warning('No nonzero thermodilution CO found in numeric file. Using uncalibrated scale.');
    k_cal = 1;
else
    TD0 = td_CO_Lmin_all(1); %first measured TCO value
    t0_sec = td_time_sec_all(1); %time of first measured CO
    [~,i_near] = min(abs(to_min - (t0_sec/60)));
    est0 = co_uncal(i_near); %uncalibrated estimate at that time
    k_cal = TD0 / est0; %Calibration factor
    fprintf('Calibration factor k = %.3f (first TD at %.2f hr)\n',k_cal,t0_sec/3600);
end
co_cal = k_cal * co_uncal; %apply calibration

% ---------------------- SAMPLE FEATURES AT TD TIMES --------------------
% extract PP, MAP, and HR values for plotting
PP_td  = []; MAP_td = []; HR_td = [];
for ii = 1:numel(td_time_sec_12h)
    [~,jj] = min(abs(to_min - (td_time_sec_12h(ii)/60)));
    PP_td(ii,1)  = PP(jj);
    MAP_td(ii,1) = MAP(jj);
    HR_td(ii,1)  = HR(jj);
end

% ---------------------- PLOT RESULTS -----------------------------------
fprintf('Plotting calibrated CO, PP, MAP, HR...\n');
figure('Color','w','Name','Problem3 Liljestrand');
sgtitle('Problem 3 – Patient #20 | Liljestrand (#5) | First 12h');


% --- 1) Continuous CO (calibrated) ---
subplot(4,1,1); hold on;
plot(to_hr,co_cal,'Color',[0.5 0.5 0.5],'LineWidth',1.3);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h,td_CO_Lmin_12h,'k','filled');
end
ylabel('CO [L/min]');
xlim([0 TMAX_HR]); grid on;
title(sprintf('Calibrated CO (k = %.3f)',k_cal));


% --- 2) Pulse Pressure (PP) ---
subplot(4,1,2); hold on;
plot(to_hr,PP,'Color',[0.5 0.5 0.5]);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h,PP_td,'k','filled');
end
ylabel('PP [mmHg]'); xlim([0 TMAX_HR]); grid on;


% --- 3) Mean Arterial Pressure (MAP) ---
subplot(4,1,3); hold on;
plot(to_hr,MAP,'Color',[0.5 0.5 0.5]);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h,MAP_td,'k','filled');
end
ylabel('MAP [mmHg]'); xlim([0 TMAX_HR]); grid on;


% --- 4) Heart Rate (HR) ---
subplot(4,1,4); hold on;
plot(to_hr,HR,'Color',[0.5 0.5 0.5]);
if ~isempty(td_time_hr_12h)
    stem(td_time_hr_12h,HR_td,'k','filled');
end
ylabel('HR [bpm]'); xlabel('Time [hours]');
xlim([0 TMAX_HR]); grid on;

% Save figure to OUT_DIR
ax = findall(gcf,'type','axes');
set(ax, 'Toolbar', []);
saveas(gcf,fullfile(OUT_DIR,'Problem3_patient20.png'));
fprintf('Saved Figure to %s\n',OUT_DIR);
fprintf('\n=== Problem 3 Completed===\n');

%% ======================================================================
%% -------------------------- PROBLEM 5 ---------------------------------
%% Continuous CO estimation using Parlikar (2007) estimator (#14)
%% Three patients, auto subdirs, Problem 4-style plots
%% ======================================================================

fprintf('\n=== Starting Problem 5: Parlikar (#14) for three patients ===\n');

% ---------------------- USER-EDIT THIS BLOCK ---------------------------
ABP_FILES = { ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s00020-2567-03-30-17-47_ABP.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10013-2564-11-01-23-37_ABP.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10205-2631-06-14-11-38_ABP.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10241-3150-11-15-13-27_ABP.txt' ...
};
NUM_FILES = { ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s00020-2567-03-30-17-47n.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10013-2564-11-01-23-37n.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10205-2631-06-14-11-38n.txt', ...
    '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/main/s10241-3150-11-15-13-27n.txt' ...
};
OUT_ROOT  = '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/Q5';
% ----------------------------------------------------------------------

Fs = 125;
TMAX_HR  = 12;
TMAX_SEC = TMAX_HR * 3600;
if ~exist(OUT_ROOT,'dir'); mkdir(OUT_ROOT); end

for p = 1:numel(ABP_FILES)
    PATH_ABP = ABP_FILES{p};
    PATH_NUM = NUM_FILES{p};

    % Per-patient subdirectory
    [~, abp_base, ~] = fileparts(PATH_ABP);
    pat_label = sprintf('patient_%02d_%s', p, abp_base);
    OUT_DIR = fullfile(OUT_ROOT, pat_label);
    if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end
    fprintf('\n================ %s =================\n', pat_label);

    % ------------------------------- READ DATA -------------------------
    raw = readmatrix(PATH_ABP,'FileType','text'); % [time(s), ABP(mmHg)]
    if size(raw,2) < 2, error('ABP file must have at least 2 columns [time, ABP].'); end
    t   = raw(:,1);   abp = raw(:,2);  t = t(:);  abp = abp(:);

    % ------------------------- TRIM ABP (0–12h) ------------------------
    mask12 = (t <= TMAX_SEC);
    time12 = t(mask12);
    ABP12  = abp(mask12);

    % ---------------------- RUN 2analyze FUNCTIONS ---------------------
    fprintf('Running WABP / ABPFEATURE / JSQL for 12-hour data (Q5)\n');
    t_on  = wabp(ABP12);
    feat  = abpfeature(ABP12, t_on);
    beatq = jSQI(feat, t_on(1:end-1), ABP12);
    if size(feat,1) - 1 ~= size(beatq,1)
        M = min(size(beatq,1), size(feat,1)-1);
        beatq = beatq(1:M,:);
    end

    % Save for estimator wrapper
    time = time12; ABP = ABP12;
    mat_out_q5 = fullfile(OUT_DIR, 'first12h_Q5.mat');
    save(mat_out_q5,'time','ABP','t_on','feat','beatq');

    % ---------------------- ESTIMATE CO (UNCALIBRATED) -----------------
    fprintf('Running Parlikar estimator (#14)\n');
    [co_uncal,to_min,~,fea] = estimateCO_v2(mat_out_q5, 14, 1);  % <== capture fea
    to_hr = to_min / 60;

    % Pull features (match Part 3/4)
    Psys   = fea(:,2);
    Pdias  = fea(:,4);
    PP     = fea(:,5);
    MAP    = fea(:,6);
    Period = fea(:,7);
    HR     = 60 * Fs ./ Period;

    % ---------------------- LOAD THERMODILUTION CO ---------------------
    fprintf('Reading numeric file for CO measurements\n');
    tbl = readtable(PATH_NUM,'FileType','text','HeaderLines',2);
    time_num_all = tbl{:,1};   % seconds
    CO_TD_all    = tbl{:,end}; % CO [L/min]

    nz_idx_all      = find(CO_TD_all ~= 0 & ~isnan(CO_TD_all));
    td_time_sec_all = time_num_all(nz_idx_all);
    td_CO_Lmin_all  = CO_TD_all(nz_idx_all);

    mask12_TD       = td_time_sec_all <= TMAX_SEC;
    td_time_sec_12h = td_time_sec_all(mask12_TD);
    td_CO_Lmin_12h  = td_CO_Lmin_all(mask12_TD);
    td_time_hr_12h  = td_time_sec_12h / 3600;

    % ---------------------- CALIBRATION (C2/C3 style) ------------------
    if isempty(td_time_sec_all)
        warning('No thermodilution CO found. Using k=1.');
        k_cal = 1;
    else
        TD0    = td_CO_Lmin_all(1);
        t0_sec = td_time_sec_all(1);
        [~, i_near] = min(abs(to_min - (t0_sec/60)));
        if isempty(i_near) || i_near<1 || i_near>numel(co_uncal)
            k_cal = 1;
            warning('Calibration: invalid nearest estimator sample; default k=1.');
        else
            est0  = co_uncal(i_near);
            k_cal = TD0 / est0;
        end
        fprintf('Calibration factor k = %.3f (first TD at %.2f hr)\n', k_cal, t0_sec/3600);
    end
    co_cal = k_cal * co_uncal;

    % ---------------------- ALIGNMENT GUARD ----------------------------
    to_hr  = to_hr(:);
    co_cal = co_cal(:);
    L = min(length(to_hr), length(co_cal));
    to_hr  = to_hr(1:L);
    co_cal = co_cal(1:L);

    % ===== sample PP/MAP/HR at TD times (for stems in lower panels) ====
    PP_td  = []; MAP_td = []; HR_td = [];
    if ~isempty(td_time_sec_12h)
        for ii = 1:numel(td_time_sec_12h)
            [~,jj] = min(abs(to_min - (td_time_sec_12h(ii)/60)));
            PP_td(ii,1)  = PP(jj);
            MAP_td(ii,1) = MAP(jj);
            HR_td(ii,1)  = HR(jj);
        end
    end

    % ---------------------- EXISTING PLOTS (Problem 4 style) -----------
    % Individual algorithm plot
    figure('Color','w','Name',sprintf('Parlikar - %s', pat_label));
    hold on;
    plot(to_hr, co_cal, 'b-', 'LineWidth', 2, 'DisplayName','Parlikar (#14)');
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2, ...
             'DisplayName','Thermodilution CO');
    end
    ylabel('CO [L/min]'); xlabel('Time [hours]');
    xlim([0 TMAX_HR]); grid on;
    title(sprintf('Parlikar Algorithm (#14) - %s', pat_label), 'Interpreter','none');
    legend('Location','best');
    saveas(gcf, fullfile(OUT_DIR, 'Problem5_Parlikar.png'));

    % Comparison-style plot
    figure('Color','w','Name',sprintf('Comparison - %s', pat_label));
    hold on;
    plot(to_hr, co_cal, 'b-', 'LineWidth', 1.5, 'DisplayName','Parlikar (#14)');
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled', 'LineWidth', 2, ...
             'DisplayName','Thermodilution CO');
    end
    ylabel('CO [L/min]'); xlabel('Time [hours]');
    xlim([0 TMAX_HR]); grid on; legend('Location','best');
    title(sprintf('Algorithm Comparison - Parlikar (#14) vs TD | %s', pat_label), 'Interpreter','none');
    saveas(gcf, fullfile(OUT_DIR, 'Problem5_Parlikar_Comparison.png'));

    % ---------------------- NEW: Figure-4 style (paper) -----------------
    % Grey lines, black stems in all 4 panels: CO, PP, MAP, HR
    grey = [0.5 0.5 0.5];

    f4 = figure('Color','w','Name',sprintf('Figure4 Style - %s', pat_label));
    tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

    % --- CO panel ---
    nexttile; hold on;
    plot(to_hr, co_cal, 'Color', grey, 'LineWidth', 1);
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, td_CO_Lmin_12h, 'k', 'filled');
    end
    ylabel('CO'); xlim([0 TMAX_HR]); grid on;

    % --- PP panel ---
    nexttile; hold on;
    plot(to_hr(1:min(L,numel(PP))), PP(1:min(L,numel(PP))), 'Color', grey, 'LineWidth', 1);
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, PP_td, 'k', 'filled');
    end
    ylabel('PP'); xlim([0 TMAX_HR]); grid on;

    % --- MAP panel ---
    nexttile; hold on;
    plot(to_hr(1:min(L,numel(MAP))), MAP(1:min(L,numel(MAP))), 'Color', grey, 'LineWidth', 1);
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, MAP_td, 'k', 'filled');
    end
    ylabel('MAP'); xlim([0 TMAX_HR]); grid on;

    % --- HR panel ---
    nexttile; hold on;
    plot(to_hr(1:min(L,numel(HR))), HR(1:min(L,numel(HR))), 'Color', grey, 'LineWidth', 1);
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h, HR_td, 'k', 'filled');
    end
    ylabel('HR'); xlabel('time [hours]');
    xlim([0 TMAX_HR]); grid on;

    % Save Figure-4 style output
    saveas(f4, fullfile(OUT_DIR, 'Problem5_Parlikar_Figure4Style.png'));
    close(f4);

    fprintf('Saved figures to %s\n', OUT_DIR);
end

compare_algorithms_q5(ABP_FILES, NUM_FILES, OUT_ROOT, [14 5 2 7]);
fprintf('\n=== Problem 5 Completed for All Three Patients ===\n');

%% -------------------------- PROBLEM 6 ---------------------------------
%% Estimate and plot TPR (Resistance) over first 12 h, Parlikar style
%   Output: one Q6 folder under OUT_ROOT with per-patient subfolders
%   Units: TPR in mmHg/(mL/s)
%   Patients: #20, 10205, 10241  (indices [1 3 4])

fprintf('\n=== Starting Problem 6: TPR Parlikar-style single-panel plots ===\n');

Fs       = 125;                % Hz
TMAX_HR  = 12;
TMAX_SEC = TMAX_HR * 3600;

% Separate Q6 root under the general output folder
OUT_Q6_ROOT = fullfile(OUT_ROOT, 'Q6');
if ~exist(OUT_Q6_ROOT,'dir'); mkdir(OUT_Q6_ROOT); end

% Patients to process
PATS_FOR_Q6 = [1 3 4];

for idx = 1:numel(PATS_FOR_Q6)
    p = PATS_FOR_Q6(idx);

    PATH_ABP = ABP_FILES{p};
    PATH_NUM = NUM_FILES{p};
    [~, abp_base, ~] = fileparts(PATH_ABP);
    pat_label = sprintf('patient_%02d_%s', p, abp_base);

    % Problem 5 data and new Q6 folder
    OUT_DIR_Q5 = fullfile(OUT_ROOT, pat_label);
    OUT_DIR_Q6 = fullfile(OUT_Q6_ROOT, pat_label);
    if ~exist(OUT_DIR_Q6,'dir'); mkdir(OUT_DIR_Q6); end

    fprintf('\n--- Generating TPR plot for %s ---\n', pat_label);

    % Load saved 12-hour ABP features
    mat_out_q5 = fullfile(OUT_DIR_Q5, 'first12h_Q5.mat');
    if ~isfile(mat_out_q5)
        error('Missing %s. Run Problem 5 first.', mat_out_q5);
    end

    % Run Parlikar estimator (#14)
    [co_uncal, to_min, ~, fea] = estimateCO_v2(mat_out_q5, 14, 1);
    to_sec = to_min * 60;

    % Extract needed features
    MAP    = fea(:,6);               % mean ABP [mmHg]
    Period = fea(:,7);
    HR     = 60 * Fs ./ Period; %#ok<NASGU>

    % ---- Calibration using first thermodilution point ----
    tbl = readtable(PATH_NUM,'FileType','text','HeaderLines',2);
    time_num_all = tbl{:,1};
    CO_TD_all    = tbl{:,end};
    nz_idx_all      = find(CO_TD_all ~= 0 & ~isnan(CO_TD_all));
    td_time_sec_all = time_num_all(nz_idx_all);
    td_CO_Lmin_all  = CO_TD_all(nz_idx_all);

    if isempty(td_time_sec_all)
        k_cal = 1;
    else
        TD0 = td_CO_Lmin_all(1);
        t0_sec = td_time_sec_all(1);
        [~, i_near] = min(abs(to_min - (t0_sec/60)));
        if isempty(i_near) || i_near<1 || i_near>numel(co_uncal)
            k_cal = 1;
        else
            k_cal = TD0 / co_uncal(i_near);
        end
    end
    co_cal = k_cal * co_uncal;

    % Align vectors
    L = min([length(to_sec), length(co_cal), length(MAP)]);
    to_sec = to_sec(1:L);
    co_cal = co_cal(1:L);
    MAP    = MAP(1:L);

    % ---- Compute TPR ----
    CO_mlps = co_cal * 1000 / 60;       % [mL/s]
    TPR = MAP ./ CO_mlps;               % [mmHg/(mL/s)]
    TPR(~isfinite(TPR)) = NaN;

    % ---- Plot Parlikar-style single-panel Resistance curve ----
    f = figure('Color','w','Name',sprintf('TPR_%s', pat_label));
    plot(to_sec, TPR, 'b', 'LineWidth', 1.2); hold on;
    ylabel('Resistance (mmHg/(mL/s))');
    xlabel('Time (secs)');
    title(sprintf('Estimated TPR – %s', pat_label), 'Interpreter','none');
    grid on; xlim([0 TMAX_SEC]);

    % Save to Q6 subdirectory
    saveas(f, fullfile(OUT_DIR_Q6, sprintf('%s_TPR_only.png', pat_label)));
    close(f);

    fprintf('Saved Parlikar-style TPR plot to %s\n', OUT_DIR_Q6);
end

fprintf('\n=== Problem 6 (TPR-only) Completed for Selected Patients ===\n');