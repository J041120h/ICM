function task4_popPD
% Task 4 – Population PD:
% Vary tumor VEGF and VEGFR2 expression using VirtualPatientsTask4.csv
% and simulate 100 virtual patients under a metronomic regimen
% (1 mg/kg daily × 10 days, total 30 days).
%
% Visualizes VEGF, VEGFR2, and VEGF–VEGFR2 in:
%   - blood
%   - tumor
%   - rest-of-body ("normal tissue")
%
% Outputs: task4_popPD.png

    clearvars -except task4_popPD;
    clc;

    %% Index mapping (same as Tasks 1–2)
    n.VEGF_b        = 1;
    n.VEGFR2_b      = 2;
    n.VEGFVEGFR2_b  = 3;
    n.Ab_b          = 4;
    n.VEGFAb_b      = 5;

    n.VEGF_t        = 6;
    n.VEGFR2_t      = 7;
    n.VEGFVEGFR2_t  = 8;
    n.Ab_t          = 9;
    n.VEGFAb_t      = 10;

    n.VEGF_r        = 11;
    n.VEGFR2_r      = 12;
    n.VEGFVEGFR2_r  = 13;
    n.Ab_r          = 14;
    n.VEGFAb_r      = 15;

    %% Baseline parameters p (same as in tasks 1–2)
    p.Vol_b = 5 * 0.60;    % L, plasma volume of blood
    p.Vol_t = 1 * 0.611;   % L, accessible tumor volume
    p.Vol_r = 40 * 0.0816; % L, accessible rest-of-body volume

    % Clearance (min^-1)
    p.kcl_V   = 0.0648;
    p.kcl_VA  = 2.2e-5;
    p.kcl_A   = 2.2e-5;

    % Internalization
    p.kintR2  = 0.0168;
    p.kintVR2 = 0.0168;

    % VEGF production (pM/min)
    p.VEGFprod_b   = 0;
    p.VEGFprod_t   = 3;   % baseline tumor VEGF production
    p.VEGFprod_r   = 10;

    % VEGFR2 production – baseline from ss conc * kintR2
    p.VEGFR2prod_b = 0;
    p.VEGFR2prod_t = 2.85e2 * p.kintR2;   % tumor
    p.VEGFR2prod_r = 2.2e3  * p.kintR2;   % rest of body

    % Binding parameters
    p.konVR  = 6e-4;
    p.koffVR = 0.06;

    p.konVA  = 5.52e-6;
    p.koffVA = 0.012;

    % Transport parameters
    p.k_bt = 4.12e-4;
    p.k_br = 3.19e-3;
    p.k_tb = 4.12e-4;
    p.k_rb = 3.19e-3 + 6.15e-4;

    % Allow antibody extravasation
    p.AbEx = 1;

    %% Dose – same as Task 2, then scaled to 1 mg/kg daily
    dose10_mg_total   = 10 * 70;    % mg
    MW_mg_per_mol     = 150000000;  % mg/mol (150 kDa)
    dose10_mol_total  = dose10_mg_total / MW_mg_per_mol;
    dose10_pmol_total = dose10_mol_total * 1e12;  % pmol
    dose10_conc_pM    = dose10_pmol_total / p.Vol_b;

    % Metronomic: 1 mg/kg daily × 10 days = 1/10 of 10 mg/kg daily
    dose1_conc_pM = dose10_conc_pM / 10;

    %% First: get steady state without antibody
    y0      = zeros(15,1);
    sstime  = 60 * 24 * 10;       % 10 days to reach steady state
    tspanSS = 0:1:sstime;         % 1-min steps

    solverSS  = @ode15s;
    optionsSS = odeset('MaxStep', 5e-1, ...
                       'AbsTol', 1e-5, ...
                       'RelTol', 1e-5, ...
                       'InitialStep', 1e-2);

    [~, Yss] = solverSS(@VEGFAbeqns, tspanSS, y0, optionsSS, p, n);
    y_ss = Yss(end, :).';

    %% Time settings for metronomic regimen
    n_days_dose     = 10;              % daily dosing for 10 days
    total_days      = 30;              % simulate out to 30 days
    minutes_per_day = 24 * 60;
    final_time      = total_days * minutes_per_day;

    %% Read virtual patient PD data from CSV
    vpTable  = readtable('VirtualPatientsTask4.csv');
    varNames = vpTable.Properties.VariableNames;

    idx_VEGF   = find(contains(varNames, 'sVEGF'),   1);
    idx_VEGFR2 = find(contains(varNames, 'sVEGFR2'), 1);

    if isempty(idx_VEGF) || isempty(idx_VEGFR2)
        error('Could not find sVEGF or sVEGFR2 columns in VirtualPatientsTask4.csv');
    end

    sVEGF_tumor   = vpTable{:, idx_VEGF};
    sVEGFR2_tumor = vpTable{:, idx_VEGFR2};

    nPatients = height(vpTable);

    %% Run first patient (outside parfor) to define common time grid
    p_first = p;
    p_first.VEGFprod_t   = sVEGF_tumor(1);              % pM/min directly
    p_first.VEGFR2prod_t = sVEGFR2_tumor(1) * p.kintR2; % ss conc * kintR2

    [T_base, Y_base] = run_metronomic_for_p( ...
        p_first, y_ss, minutes_per_day, final_time, n, dose1_conc_pM, n_days_dose);

    nTime = numel(T_base);

    % Preallocate for all compartments
    VEGF_b_all        = zeros(nTime, nPatients);
    VEGF_t_all        = zeros(nTime, nPatients);
    VEGF_r_all        = zeros(nTime, nPatients);

    VEGFR2_b_all      = zeros(nTime, nPatients);
    VEGFR2_t_all      = zeros(nTime, nPatients);
    VEGFR2_r_all      = zeros(nTime, nPatients);

    VEGFVEGFR2_b_all  = zeros(nTime, nPatients);
    VEGFVEGFR2_t_all  = zeros(nTime, nPatients);
    VEGFVEGFR2_r_all  = zeros(nTime, nPatients);

    % Fill first patient
    VEGF_b_all(:,1)       = Y_base(:, n.VEGF_b);
    VEGF_t_all(:,1)       = Y_base(:, n.VEGF_t);
    VEGF_r_all(:,1)       = Y_base(:, n.VEGF_r);

    VEGFR2_b_all(:,1)     = Y_base(:, n.VEGFR2_b);
    VEGFR2_t_all(:,1)     = Y_base(:, n.VEGFR2_t);
    VEGFR2_r_all(:,1)     = Y_base(:, n.VEGFR2_r);

    VEGFVEGFR2_b_all(:,1) = Y_base(:, n.VEGFVEGFR2_b);
    VEGFVEGFR2_t_all(:,1) = Y_base(:, n.VEGFVEGFR2_t);
    VEGFVEGFR2_r_all(:,1) = Y_base(:, n.VEGFVEGFR2_r);

    %% Run remaining patients in parallel
    % If parfor causes issues, change to `for`.
    parfor i = 2:nPatients
        p_i = p;
        p_i.VEGFprod_t   = sVEGF_tumor(i);
        p_i.VEGFR2prod_t = sVEGFR2_tumor(i) * p.kintR2;

        [T_i, Y_i] = run_metronomic_for_p( ...
            p_i, y_ss, minutes_per_day, final_time, n, dose1_conc_pM, n_days_dose);

        VEGF_b_all(:,i)       = Y_i(:, n.VEGF_b);
        VEGF_t_all(:,i)       = Y_i(:, n.VEGF_t);
        VEGF_r_all(:,i)       = Y_i(:, n.VEGF_r);

        VEGFR2_b_all(:,i)     = Y_i(:, n.VEGFR2_b);
        VEGFR2_t_all(:,i)     = Y_i(:, n.VEGFR2_t);
        VEGFR2_r_all(:,i)     = Y_i(:, n.VEGFR2_r);

        VEGFVEGFR2_b_all(:,i) = Y_i(:, n.VEGFVEGFR2_b);
        VEGFVEGFR2_t_all(:,i) = Y_i(:, n.VEGFVEGFR2_t);
        VEGFVEGFR2_r_all(:,i) = Y_i(:, n.VEGFVEGFR2_r);
    end

    %% Plot VEGF / VEGFR2 / VEGF–VEGFR2 in blood, tumor, rest-of-body
    t_days = T_base / (60*24);
    lw_ind = 0.4;
    lw_mean = 2;

    fig = figure('Visible','off');
    fig.Position = [100 100 1400 900];

    % Row 1: VEGF
    subplot(3,3,1);
    hold on;
    plot(t_days, VEGF_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_b_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Blood');
    xlabel('Time (days)');
    ylabel('[VEGF]_b (pM)');

    subplot(3,3,2);
    hold on;
    plot(t_days, VEGF_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_t_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Tumor');
    xlabel('Time (days)');
    ylabel('[VEGF]_t (pM)');

    subplot(3,3,3);
    hold on;
    plot(t_days, VEGF_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_r_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Rest of body');
    xlabel('Time (days)');
    ylabel('[VEGF]_r (pM)');

    % Row 2: VEGFR2
    subplot(3,3,4);
    hold on;
    plot(t_days, VEGFR2_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFR2_b_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGFR2 – Blood');
    xlabel('Time (days)');
    ylabel('[VEGFR2]_b (pM)');

    subplot(3,3,5);
    hold on;
    plot(t_days, VEGFR2_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFR2_t_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGFR2 – Tumor');
    xlabel('Time (days)');
    ylabel('[VEGFR2]_t (pM)');

    subplot(3,3,6);
    hold on;
    plot(t_days, VEGFR2_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFR2_r_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGFR2 – Rest of body');
    xlabel('Time (days)');
    ylabel('[VEGFR2]_r (pM)');

    % Row 3: VEGF–VEGFR2 complex
    subplot(3,3,7);
    hold on;
    plot(t_days, VEGFVEGFR2_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFVEGFR2_b_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–VEGFR2 – Blood');
    xlabel('Time (days)');
    ylabel('[VEGF–VEGFR2]_b (pM)');

    subplot(3,3,8);
    hold on;
    plot(t_days, VEGFVEGFR2_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFVEGFR2_t_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–VEGFR2 – Tumor');
    xlabel('Time (days)');
    ylabel('[VEGF–VEGFR2]_t (pM)');

    subplot(3,3,9);
    hold on;
    plot(t_days, VEGFVEGFR2_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFVEGFR2_r_all,2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–VEGFR2 – Rest of body');
    xlabel('Time (days)');
    ylabel('[VEGF–VEGFR2]_r (pM)');

    sgtitle('Task 4 – PopPD: VEGF/VEGFR2 Variability in Blood, Tumor, and Normal Tissue');

    saveas(fig, 'task4_popPD.png');
    close(fig);
end

% =====================================================================
% Subfunction: run metronomic regimen for a given parameter set p_local
% =====================================================================
function [T_all, Y_all] = run_metronomic_for_p( ...
    p_local, y_ss, minutes_per_day, final_time, n, dose1_conc_pM, n_days_dose)

    solver  = @ode15s;
    options = odeset('MaxStep', 5e-1, ...
                     'AbsTol', 1e-5, ...
                     'RelTol', 1e-5, ...
                     'InitialStep', 1e-2);

    y_curr = y_ss;
    t_curr = 0;

    T_all = [];
    Y_all = [];

    % Daily dosing phase
    for d = 1:n_days_dose
        % Give that day's dose into blood
        y_curr(n.Ab_b) = y_curr(n.Ab_b) + dose1_conc_pM;

        tspan_day = t_curr:1:(t_curr + minutes_per_day);
        [T_day, Y_day] = solver(@VEGFAbeqns, tspan_day, y_curr, options, p_local, n);

        if isempty(T_all)
            T_all = T_day;
            Y_all = Y_day;
        else
            T_all = [T_all; T_day(2:end)];
            Y_all = [Y_all; Y_day(2:end,:)];
        end

        t_curr = t_curr + minutes_per_day;
        y_curr = Y_day(end, :).';
    end

    % Tail phase: no further doses, simulate until final_time
    if t_curr < final_time
        tspan_tail = t_curr:1:final_time;
        [T_tail, Y_tail] = solver(@VEGFAbeqns, tspan_tail, y_curr, options, p_local, n);
        T_all = [T_all; T_tail(2:end)];
        Y_all = [Y_all; Y_tail(2:end,:)];
    end
end
