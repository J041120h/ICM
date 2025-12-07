function task4_popPD
% Task 4 – Population PD:
% Vary tumor VEGF and VEGFR2 expression using VirtualPatientsTask4.csv
% and simulate 100 virtual patients under:
%   - Single injection (10 mg/kg, once at t = 0)
%   - Metronomic regimen (1 mg/kg daily × 10 days)
% For each regimen, run two scenarios:
%   - No antibody extravasation (AbEx = 0)
%   - With antibody extravasation (AbEx = 1)
%
% For each of the 4 scenarios, generate one 3x5 multipanel figure:
%   Rows:   blood, tumor, rest of body
%   Cols:   VEGF, VEGFR2, VEGF–VEGFR2, Ab, VEGF–Ab
%
% Total: 4 figures × (3*5 panels) = 60 subplots.

    clearvars -except task4_popPD;
    clc;

    %% Index mapping (consistent with VEGFAbeqns.m)
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

    %% Baseline parameters p (same as Tasks 1–2, except AbEx will vary by scenario)
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

    % We will overwrite p.AbEx for each scenario (0 or 1)
    p.AbEx = 1;  % placeholder; irrelevant for steady state (Ab=0)

    %% Dose calculations
    % Single-injection: 10 mg/kg (70 kg patient)
    dose10_mg_total   = 10 * 70;          % mg
    MW_mg_per_mol     = 150000000;        % mg/mol (150 kDa)
    dose10_mol_total  = dose10_mg_total / MW_mg_per_mol;
    dose10_pmol_total = dose10_mol_total * 1e12;  % pmol
    dose10_conc_pM    = dose10_pmol_total / p.Vol_b;

    % Metronomic: 1 mg/kg daily × 10 days = 1/10 of 10 mg/kg daily
    dose1_conc_pM = dose10_conc_pM / 10;

    %% First: get steady state without antibody (Ab = 0 everywhere)
    y0      = zeros(15,1);
    sstime  = 60 * 24 * 10;       % 10 days to reach steady state
    tspanSS = 0:1:sstime;         % 1-min steps

    solverSS  = @ode15s;
    optionsSS = odeset('MaxStep',    5e-1, ...
                       'AbsTol',     1e-5, ...
                       'RelTol',     1e-5, ...
                       'InitialStep',1e-2);

    [~, Yss] = solverSS(@VEGFAbeqns, tspanSS, y0, optionsSS, p, n);
    y_ss = Yss(end, :).';

    %% Time settings (same for single and metronomic)
    n_days_dose_metro = 10;              % daily dosing for 10 days (metronomic only)
    total_days        = 30;              % simulate out to 30 days
    minutes_per_day   = 24 * 60;
    final_time        = total_days * minutes_per_day;

    %% Read virtual patient PD data from CSV (tumor VEGF and VEGFR2 production)
    vpTable  = readtable('VirtualPatientsTask4.csv');  % relative path
    varNames = vpTable.Properties.VariableNames;

    idx_VEGF   = find(contains(varNames, 'sVEGF'),   1);
    idx_VEGFR2 = find(contains(varNames, 'sVEGFR2'), 1);

    if isempty(idx_VEGF) || isempty(idx_VEGFR2)
        error('Could not find sVEGF or sVEGFR2 columns in VirtualPatientsTask4.csv');
    end

    sVEGF_tumor   = vpTable{:, idx_VEGF};
    sVEGFR2_tumor = vpTable{:, idx_VEGFR2};

    nPatients = height(vpTable);

    %% Define the 4 scenarios: regimen × AbEx
    scenarios = struct( ...
        'label',   {}, ...
        'title',   {}, ...
        'regimen', {}, ...
        'AbEx',    [] , ...
        'outfile', {}  );

    scenarios(1).label   = 'single_noAbEx';
    scenarios(1).title   = 'Single injection (10 mg/kg) – No Ab extravasation';
    scenarios(1).regimen = 'single';
    scenarios(1).AbEx    = 0;
    scenarios(1).outfile = 'task4_popPD_single_noAbEx.png';

    scenarios(2).label   = 'single_withAbEx';
    scenarios(2).title   = 'Single injection (10 mg/kg) – With Ab extravasation';
    scenarios(2).regimen = 'single';
    scenarios(2).AbEx    = 1;
    scenarios(2).outfile = 'task4_popPD_single_withAbEx.png';

    scenarios(3).label   = 'metronomic_noAbEx';
    scenarios(3).title   = 'Metronomic (1 mg/kg × 10 days) – No Ab extravasation';
    scenarios(3).regimen = 'metronomic';
    scenarios(3).AbEx    = 0;
    scenarios(3).outfile = 'task4_popPD_metronomic_noAbEx.png';

    scenarios(4).label   = 'metronomic_withAbEx';
    scenarios(4).title   = 'Metronomic (1 mg/kg × 10 days) – With Ab extravasation';
    scenarios(4).regimen = 'metronomic';
    scenarios(4).AbEx    = 1;
    scenarios(4).outfile = 'task4_popPD_metronomic_withAbEx.png';

    %% Line style settings
    lw_ind  = 0.4;   % individual patients
    lw_mean = 2.0;   % population mean

    %% Loop over each scenario and generate one 3x5 figure
    for sIdx = 1:numel(scenarios)
        scn = scenarios(sIdx);
        fprintf('Running scenario %d/%d: %s\n', sIdx, numel(scenarios), scn.label);

        % Scenario-specific parameters
        p_scn        = p;
        p_scn.AbEx   = scn.AbEx;

        % --- First patient (to define common time grid) ---
        p_first = p_scn;
        p_first.VEGFprod_t   = sVEGF_tumor(1);              % pM/min directly
        p_first.VEGFR2prod_t = sVEGFR2_tumor(1) * p.kintR2; % ss conc * kintR2

        [T_base, Y_base] = run_regimen_for_p( ...
            p_first, y_ss, minutes_per_day, final_time, n, ...
            dose10_conc_pM, dose1_conc_pM, n_days_dose_metro, scn.regimen);

        nTime  = numel(T_base);
        t_days = T_base / (60*24);

        % --- Preallocate: all 5 molecules × 3 compartments ---
        VEGF_b_all       = zeros(nTime, nPatients);
        VEGF_t_all       = zeros(nTime, nPatients);
        VEGF_r_all       = zeros(nTime, nPatients);

        VEGFR2_b_all     = zeros(nTime, nPatients);
        VEGFR2_t_all     = zeros(nTime, nPatients);
        VEGFR2_r_all     = zeros(nTime, nPatients);

        VEGFVEGFR2_b_all = zeros(nTime, nPatients);
        VEGFVEGFR2_t_all = zeros(nTime, nPatients);
        VEGFVEGFR2_r_all = zeros(nTime, nPatients);

        Ab_b_all         = zeros(nTime, nPatients);
        Ab_t_all         = zeros(nTime, nPatients);
        Ab_r_all         = zeros(nTime, nPatients);

        VEGFAb_b_all     = zeros(nTime, nPatients);
        VEGFAb_t_all     = zeros(nTime, nPatients);
        VEGFAb_r_all     = zeros(nTime, nPatients);

        % Fill arrays for first patient
        VEGF_b_all(:,1)       = Y_base(:, n.VEGF_b);
        VEGF_t_all(:,1)       = Y_base(:, n.VEGF_t);
        VEGF_r_all(:,1)       = Y_base(:, n.VEGF_r);

        VEGFR2_b_all(:,1)     = Y_base(:, n.VEGFR2_b);
        VEGFR2_t_all(:,1)     = Y_base(:, n.VEGFR2_t);
        VEGFR2_r_all(:,1)     = Y_base(:, n.VEGFR2_r);

        VEGFVEGFR2_b_all(:,1) = Y_base(:, n.VEGFVEGFR2_b);
        VEGFVEGFR2_t_all(:,1) = Y_base(:, n.VEGFVEGFR2_t);
        VEGFVEGFR2_r_all(:,1) = Y_base(:, n.VEGFVEGFR2_r);

        Ab_b_all(:,1)         = Y_base(:, n.Ab_b);
        Ab_t_all(:,1)         = Y_base(:, n.Ab_t);
        Ab_r_all(:,1)         = Y_base(:, n.Ab_r);

        VEGFAb_b_all(:,1)     = Y_base(:, n.VEGFAb_b);
        VEGFAb_t_all(:,1)     = Y_base(:, n.VEGFAb_t);
        VEGFAb_r_all(:,1)     = Y_base(:, n.VEGFAb_r);

        % --- Remaining patients in parallel ---
        parfor i = 2:nPatients
            p_i = p_scn;
            p_i.VEGFprod_t   = sVEGF_tumor(i);
            p_i.VEGFR2prod_t = sVEGFR2_tumor(i) * p.kintR2;

            [T_i, Y_i] = run_regimen_for_p( ...
                p_i, y_ss, minutes_per_day, final_time, n, ...
                dose10_conc_pM, dose1_conc_pM, n_days_dose_metro, scn.regimen);

            % assume T_i matches T_base (same regimen and time grid)
            VEGF_b_all(:,i)       = Y_i(:, n.VEGF_b);
            VEGF_t_all(:,i)       = Y_i(:, n.VEGF_t);
            VEGF_r_all(:,i)       = Y_i(:, n.VEGF_r);

            VEGFR2_b_all(:,i)     = Y_i(:, n.VEGFR2_b);
            VEGFR2_t_all(:,i)     = Y_i(:, n.VEGFR2_t);
            VEGFR2_r_all(:,i)     = Y_i(:, n.VEGFR2_r);

            VEGFVEGFR2_b_all(:,i) = Y_i(:, n.VEGFVEGFR2_b);
            VEGFVEGFR2_t_all(:,i) = Y_i(:, n.VEGFVEGFR2_t);
            VEGFVEGFR2_r_all(:,i) = Y_i(:, n.VEGFVEGFR2_r);

            Ab_b_all(:,i)         = Y_i(:, n.Ab_b);
            Ab_t_all(:,i)         = Y_i(:, n.Ab_t);
            Ab_r_all(:,i)         = Y_i(:, n.Ab_r);

            VEGFAb_b_all(:,i)     = Y_i(:, n.VEGFAb_b);
            VEGFAb_t_all(:,i)     = Y_i(:, n.VEGFAb_t);
            VEGFAb_r_all(:,i)     = Y_i(:, n.VEGFAb_r);
        end

        %% === Plot: 3 (compartments) × 5 (molecules) multipanel ===
        fig = figure('Visible','off');
        fig.Position = [100 100 1600 900];

        % Row labels for compartments
        rowNames = {'Blood', 'Tumor', 'Rest of body'};

        % ----- Row 1: Blood -----
        subplot(3,5,1);
        plot_panel(t_days, VEGF_b_all,   lw_ind, lw_mean);
        title('VEGF – Blood');
        ylabel('[VEGF]_b (pM)');
        xlabel('Time (days)');

        subplot(3,5,2);
        plot_panel(t_days, VEGFR2_b_all, lw_ind, lw_mean);
        title('VEGFR2 – Blood');
        ylabel('[VEGFR2]_b (pM)');
        xlabel('Time (days)');

        subplot(3,5,3);
        plot_panel(t_days, VEGFVEGFR2_b_all, lw_ind, lw_mean);
        title('VEGF–VEGFR2 – Blood');
        ylabel('[VEGF–VEGFR2]_b (pM)');
        xlabel('Time (days)');

        subplot(3,5,4);
        plot_panel(t_days, Ab_b_all, lw_ind, lw_mean);
        title('Ab – Blood');
        ylabel('[Ab]_b (pM)');
        xlabel('Time (days)');

        subplot(3,5,5);
        plot_panel(t_days, VEGFAb_b_all, lw_ind, lw_mean);
        title('VEGF–Ab – Blood');
        ylabel('[VEGF–Ab]_b (pM)');
        xlabel('Time (days)');

        % ----- Row 2: Tumor -----
        subplot(3,5,6);
        plot_panel(t_days, VEGF_t_all, lw_ind, lw_mean);
        title('VEGF – Tumor');
        ylabel('[VEGF]_t (pM)');
        xlabel('Time (days)');

        subplot(3,5,7);
        plot_panel(t_days, VEGFR2_t_all, lw_ind, lw_mean);
        title('VEGFR2 – Tumor');
        ylabel('[VEGFR2]_t (pM)');
        xlabel('Time (days)');

        subplot(3,5,8);
        plot_panel(t_days, VEGFVEGFR2_t_all, lw_ind, lw_mean);
        title('VEGF–VEGFR2 – Tumor');
        ylabel('[VEGF–VEGFR2]_t (pM)');
        xlabel('Time (days)');

        subplot(3,5,9);
        plot_panel(t_days, Ab_t_all, lw_ind, lw_mean);
        title('Ab – Tumor');
        ylabel('[Ab]_t (pM)');
        xlabel('Time (days)');

        subplot(3,5,10);
        plot_panel(t_days, VEGFAb_t_all, lw_ind, lw_mean);
        title('VEGF–Ab – Tumor');
        ylabel('[VEGF–Ab]_t (pM)');
        xlabel('Time (days)');

        % ----- Row 3: Rest of body -----
        subplot(3,5,11);
        plot_panel(t_days, VEGF_r_all, lw_ind, lw_mean);
        title('VEGF – Rest of body');
        ylabel('[VEGF]_r (pM)');
        xlabel('Time (days)');

        subplot(3,5,12);
        plot_panel(t_days, VEGFR2_r_all, lw_ind, lw_mean);
        title('VEGFR2 – Rest of body');
        ylabel('[VEGFR2]_r (pM)');
        xlabel('Time (days)');

        subplot(3,5,13);
        plot_panel(t_days, VEGFVEGFR2_r_all, lw_ind, lw_mean);
        title('VEGF–VEGFR2 – Rest of body');
        ylabel('[VEGF–VEGFR2]_r (pM)');
        xlabel('Time (days)');

        subplot(3,5,14);
        plot_panel(t_days, Ab_r_all, lw_ind, lw_mean);
        title('Ab – Rest of body');
        ylabel('[Ab]_r (pM)');
        xlabel('Time (days)');

        subplot(3,5,15);
        plot_panel(t_days, VEGFAb_r_all, lw_ind, lw_mean);
        title('VEGF–Ab – Rest of body');
        ylabel('[VEGF–Ab]_r (pM)');
        xlabel('Time (days)');

        sgtitle(sprintf('Task 4 – PopPD: %s', scn.title));

        % Save figure for this scenario
        saveas(fig, scn.outfile);
        close(fig);
    end

    fprintf('Task 4 PopPD: Completed all 4 scenarios.\n');
end

% =====================================================================
% Subfunction: run one regimen for a given parameter set p_local
% =====================================================================
function [T_all, Y_all] = run_regimen_for_p( ...
    p_local, y_ss, minutes_per_day, final_time, n, ...
    dose10_conc_pM, dose1_conc_pM, n_days_metro, regimenType)

    solver  = @ode15s;
    options = odeset('MaxStep',    5e-1, ...
                     'AbsTol',     1e-5, ...
                     'RelTol',     1e-5, ...
                     'InitialStep',1e-2);

    y_curr = y_ss;
    t_curr = 0;

    T_all = [];
    Y_all = [];

    switch lower(regimenType)
        case 'single'
            % Single 10 mg/kg bolus at t = 0
            y_curr(n.Ab_b) = y_curr(n.Ab_b) + dose10_conc_pM;

            tspan = t_curr:1:final_time;
            [T_all, Y_all] = solver(@VEGFAbeqns, tspan, y_curr, options, p_local, n);

        case 'metronomic'
            % Daily 1 mg/kg bolus for n_days_metro, then tail with no further doses
            for d = 1:n_days_metro
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

        otherwise
            error('Unknown regimenType: %s', regimenType);
    end
end

% =====================================================================
% Subfunction: plotting helper for each panel
% =====================================================================
function plot_panel(t_days, data_all, lw_ind, lw_mean)
    hold on;
    plot(t_days, data_all, 'LineWidth', lw_ind);
    plot(t_days, mean(data_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    box on;
end
