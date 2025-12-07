function task3_population_PK
% Task 3 – Population PK with 4 regimens/transport scenarios
%
% Scenarios (each produces a 5x3 population figure):
%   1) Single injection, NO Ab extravasation
%   2) Single injection, WITH Ab extravasation
%   3) Metronomic dosing, WITH Ab extravasation
%   4) Metronomic dosing, NO Ab extravasation
%
% For each scenario, we simulate 100 virtual patients with variability in:
%   - kcl_A (Ab clearance), CV = 50%
%   - kcl_VA (VEGF-Ab clearance), CV = 50%
%   - Vol_r (rest-of-body volume), CV = 25%
%
% Each figure shows 5 molecule types x 3 compartments:
%   Rows (species): VEGF (pM), Ab (μM), VEGF-Ab (nM), VEGFR2 (pM), VEGF-VEGFR2 (pM)
%   Cols (compartments): Blood, Tumor, Rest-of-body

    clearvars -except task3_population_PK;
    clc;

    %% Index mapping for state vector y
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

    %% Baseline parameters
    p.Vol_b = 5 * 0.60;
    p.Vol_t = 1 * 0.611;
    p.Vol_r = 40 * 0.0816;

    p.kcl_V   = 0.0648;
    p.kcl_VA  = 2.2e-5;
    p.kcl_A   = 2.2e-5;

    p.kintR2  = 0.0168;
    p.kintVR2 = 0.0168;

    p.VEGFprod_b   = 0;
    p.VEGFprod_t   = 3;
    p.VEGFprod_r   = 10;

    p.VEGFR2prod_b = 0;
    p.VEGFR2prod_t = 2.85e2 * p.kintR2;
    p.VEGFR2prod_r = 2.2e3  * p.kintR2;

    p.konVR  = 6e-4;
    p.koffVR = 0.06;

    p.konVA  = 5.52e-6;
    p.koffVA = 0.012;

    p.k_bt = 4.12e-4;
    p.k_br = 3.19e-3;
    p.k_tb = 4.12e-4;
    p.k_rb = 3.19e-3 + 6.15e-4;

    % Default AbEx; will overwrite per scenario
    p.AbEx = 1;

    %% Dose: 10 mg/kg IV (for single injection and basis for metronomic)
    dose_mg_total   = 10 * 70;
    MW_mg_per_mol   = 150000000;
    dose_mol_total  = dose_mg_total / MW_mg_per_mol;
    dose_pmol_total = dose_mol_total * 1e12;
    dose10_conc_pM  = dose_pmol_total / p.Vol_b;   % 10 mg/kg bolus

    % Metronomic: 1 mg/kg per day x 10 days (i.e., 1/10 of 10 mg/kg)
    dose1_conc_pM   = dose10_conc_pM / 10;

    %% Solver settings (prof's speed suggestions)
    solver  = @ode15s;
    options = odeset('MaxStep',5e-1,'AbsTol',1e-5,'RelTol',1e-5,'InitialStep',1e-2);

    %% Population definition (same across all scenarios)
    N_patients = 100;
    rng(1);

    base_kcl_A  = p.kcl_A;
    base_kcl_VA = p.kcl_VA;
    base_Vol_r  = p.Vol_r;

    kcl_A_samples  = sample_positive_normal(base_kcl_A,  0.5*base_kcl_A,  N_patients);
    kcl_VA_samples = sample_positive_normal(base_kcl_VA, 0.5*base_kcl_VA, N_patients);
    Vol_r_samples  = sample_positive_normal(base_Vol_r,  0.25*base_Vol_r, N_patients);

    %% Times for steady state and regimens
    % Steady state (no Ab)
    tspanSS      = 0:1:(60*24*10);  % 10 days

    % Single injection: simulate 21 days after dose
    single_days  = 21;
    tspan_single = 0:1:(single_days*24*60);

    % Metronomic: daily dosing for 10 days, simulate total 30 days
    metro_total_days  = 30;
    metro_dose_days   = 10;
    minutes_per_day   = 24*60;

    %% Get steady-state once (no Ab, AbEx irrelevant since Ab=0)
    y0_SS = zeros(15,1);
    [~,Yss] = solver(@VEGFAbeqns, tspanSS, y0_SS, options, p, n);
    y_ss   = Yss(end, :).';

    %% Define scenarios
    scenarios = struct([]);

    scenarios(1).name      = 'Single injection – No Ab extravasation';
    scenarios(1).filename  = 'task3_single_noAbEx.png';
    scenarios(1).regimen   = 'single';
    scenarios(1).AbEx      = 0;

    scenarios(2).name      = 'Single injection – With Ab extravasation';
    scenarios(2).filename  = 'task3_single_withAbEx.png';
    scenarios(2).regimen   = 'single';
    scenarios(2).AbEx      = 1;

    scenarios(3).name      = 'Metronomic Treatment – With Ab extravasation';
    scenarios(3).filename  = 'task3_metro_withAbEx.png';
    scenarios(3).regimen   = 'metronomic';
    scenarios(3).AbEx      = 1;

    scenarios(4).name      = 'Metronomic Treatment – N Ab extravasation';
    scenarios(4).filename  = 'task3_metro_noAbEx.png';
    scenarios(4).regimen   = 'metronomic';
    scenarios(4).AbEx      = 0;

    %% Loop over scenarios
    for s = 1:numel(scenarios)
        disp(['Running scenario: ', scenarios(s).name]);

        % Common parameter template for this scenario
        p_s = p;
        p_s.AbEx = scenarios(s).AbEx;

        %-------------------------------
        % Run first patient (ip=1) to define time grid & sizes
        %-------------------------------
        p_first = p_s;
        p_first.kcl_A  = kcl_A_samples(1);
        p_first.kcl_VA = kcl_VA_samples(1);
        p_first.Vol_r  = Vol_r_samples(1);

        [T_base, Y_base] = simulate_regimen( ...
            scenarios(s).regimen, ...
            p_first, y_ss, ...
            solver, options, ...
            dose10_conc_pM, dose1_conc_pM, ...
            tspan_single, ...
            metro_dose_days, metro_total_days, minutes_per_day, n);

        nTime  = numel(T_base);
        t_days = T_base/(60*24);

        % Preallocate for this scenario
        VEGF_b_all        = zeros(nTime,N_patients);
        VEGF_t_all        = zeros(nTime,N_patients);
        VEGF_r_all        = zeros(nTime,N_patients);

        VEGFR2_b_all      = zeros(nTime,N_patients);
        VEGFR2_t_all      = zeros(nTime,N_patients);
        VEGFR2_r_all      = zeros(nTime,N_patients);

        VEGFVEGFR2_b_all  = zeros(nTime,N_patients);
        VEGFVEGFR2_t_all  = zeros(nTime,N_patients);
        VEGFVEGFR2_r_all  = zeros(nTime,N_patients);

        Ab_b_all          = zeros(nTime,N_patients);
        Ab_t_all          = zeros(nTime,N_patients);
        Ab_r_all          = zeros(nTime,N_patients);

        VEGFAb_b_all      = zeros(nTime,N_patients);
        VEGFAb_t_all      = zeros(nTime,N_patients);
        VEGFAb_r_all      = zeros(nTime,N_patients);

        % Fill in first patient
        VEGF_b_all(:,1)       = Y_base(:,n.VEGF_b);
        VEGF_t_all(:,1)       = Y_base(:,n.VEGF_t);
        VEGF_r_all(:,1)       = Y_base(:,n.VEGF_r);

        VEGFR2_b_all(:,1)     = Y_base(:,n.VEGFR2_b);
        VEGFR2_t_all(:,1)     = Y_base(:,n.VEGFR2_t);
        VEGFR2_r_all(:,1)     = Y_base(:,n.VEGFR2_r);

        VEGFVEGFR2_b_all(:,1) = Y_base(:,n.VEGFVEGFR2_b);
        VEGFVEGFR2_t_all(:,1) = Y_base(:,n.VEGFVEGFR2_t);
        VEGFVEGFR2_r_all(:,1) = Y_base(:,n.VEGFVEGFR2_r);

        Ab_b_all(:,1)         = Y_base(:,n.Ab_b);
        Ab_t_all(:,1)         = Y_base(:,n.Ab_t);
        Ab_r_all(:,1)         = Y_base(:,n.Ab_r);

        VEGFAb_b_all(:,1)     = Y_base(:,n.VEGFAb_b);
        VEGFAb_t_all(:,1)     = Y_base(:,n.VEGFAb_t);
        VEGFAb_r_all(:,1)     = Y_base(:,n.VEGFAb_r);

        %-------------------------------
        % Remaining patients in parallel
        %-------------------------------
        parfor ip = 2:N_patients
            p_i = p_s;
            p_i.kcl_A  = kcl_A_samples(ip);
            p_i.kcl_VA = kcl_VA_samples(ip);
            p_i.Vol_r  = Vol_r_samples(ip);

            [T_i, Y_i] = simulate_regimen( ...
                scenarios(s).regimen, ...
                p_i, y_ss, ...
                solver, options, ...
                dose10_conc_pM, dose1_conc_pM, ...
                tspan_single, ...
                metro_dose_days, metro_total_days, minutes_per_day, n);

            % Assume T_i == T_base in grid
            VEGF_b_all(:,ip)       = Y_i(:,n.VEGF_b);
            VEGF_t_all(:,ip)       = Y_i(:,n.VEGF_t);
            VEGF_r_all(:,ip)       = Y_i(:,n.VEGF_r);

            VEGFR2_b_all(:,ip)     = Y_i(:,n.VEGFR2_b);
            VEGFR2_t_all(:,ip)     = Y_i(:,n.VEGFR2_t);
            VEGFR2_r_all(:,ip)     = Y_i(:,n.VEGFR2_r);

            VEGFVEGFR2_b_all(:,ip) = Y_i(:,n.VEGFVEGFR2_b);
            VEGFVEGFR2_t_all(:,ip) = Y_i(:,n.VEGFVEGFR2_t);
            VEGFVEGFR2_r_all(:,ip) = Y_i(:,n.VEGFVEGFR2_r);

            Ab_b_all(:,ip)         = Y_i(:,n.Ab_b);
            Ab_t_all(:,ip)         = Y_i(:,n.Ab_t);
            Ab_r_all(:,ip)         = Y_i(:,n.Ab_r);

            VEGFAb_b_all(:,ip)     = Y_i(:,n.VEGFAb_b);
            VEGFAb_t_all(:,ip)     = Y_i(:,n.VEGFAb_t);
            VEGFAb_r_all(:,ip)     = Y_i(:,n.VEGFAb_r);
        end

        %-------------------------------
        % Build 5x3 figure for this scenario
        %-------------------------------
        fig = figure('Visible','off');
        set(fig,'Position',[100 50 1400 1000]);

        lw_ind  = 0.4;
        lw_mean = 2;

        % ---- Row 1: VEGF (pM, y-lim [0 220]) ----
        subplot(5,3,1);
        hold on;
        plot(t_days, VEGF_b_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGF_b_all,2),'k','LineWidth',lw_mean);
        ylim([0 220]);
        title('VEGF – Blood');
        ylabel('[VEGF]_b (pM)');
        xlabel('Time (days)');

        subplot(5,3,2);
        hold on;
        plot(t_days, VEGF_t_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGF_t_all,2),'k','LineWidth',lw_mean);
        ylim([0 220]);
        title('VEGF – Tumor');
        ylabel('[VEGF]_t (pM)');
        xlabel('Time (days)');

        subplot(5,3,3);
        hold on;
        plot(t_days, VEGF_r_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGF_r_all,2),'k','LineWidth',lw_mean);
        ylim([0 220]);
        title('VEGF – Rest of Body');
        ylabel('[VEGF]_r (pM)');
        xlabel('Time (days)');

        % ---- Row 2: Ab (μM) ----
        subplot(5,3,4);
        hold on;
        plot(t_days, Ab_b_all/1e6,'LineWidth',lw_ind);
        plot(t_days, mean(Ab_b_all,2)/1e6,'k','LineWidth',lw_mean);
        title('Ab – Blood');
        ylabel('[Ab]_b (\muM)');
        xlabel('Time (days)');

        subplot(5,3,5);
        hold on;
        plot(t_days, Ab_t_all/1e6,'LineWidth',lw_ind);
        plot(t_days, mean(Ab_t_all,2)/1e6,'k','LineWidth',lw_mean);
        title('Ab – Tumor');
        ylabel('[Ab]_t (\muM)');
        xlabel('Time (days)');

        subplot(5,3,6);
        hold on;
        plot(t_days, Ab_r_all/1e6,'LineWidth',lw_ind);
        plot(t_days, mean(Ab_r_all,2)/1e6,'k','LineWidth',lw_mean);
        title('Ab – Rest of Body');
        ylabel('[Ab]_r (\muM)');
        xlabel('Time (days)');

        % ---- Row 3: VEGF–Ab (nM) ----
        subplot(5,3,7);
        hold on;
        plot(t_days, VEGFAb_b_all/1000,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFAb_b_all,2)/1000,'k','LineWidth',lw_mean);
        title('VEGF–Ab – Blood');
        ylabel('[VEGF–Ab]_b (nM)');
        xlabel('Time (days)');

        subplot(5,3,8);
        hold on;
        plot(t_days, VEGFAb_t_all/1000,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFAb_t_all,2)/1000,'k','LineWidth',lw_mean);
        title('VEGF–Ab – Tumor');
        ylabel('[VEGF–Ab]_t (nM)');
        xlabel('Time (days)');

        subplot(5,3,9);
        hold on;
        plot(t_days, VEGFAb_r_all/1000,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFAb_r_all,2)/1000,'k','LineWidth',lw_mean);
        title('VEGF–Ab – Rest of Body');
        ylabel('[VEGF–Ab]_r (nM)');
        xlabel('Time (days)');

        % ---- Row 4: VEGFR2 (pM) ----
        subplot(5,3,10);
        hold on;
        plot(t_days, VEGFR2_b_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFR2_b_all,2),'k','LineWidth',lw_mean);
        title('VEGFR2 – Blood');
        ylabel('[VEGFR2]_b (pM)');
        xlabel('Time (days)');

        subplot(5,3,11);
        hold on;
        plot(t_days, VEGFR2_t_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFR2_t_all,2),'k','LineWidth',lw_mean);
        title('VEGFR2 – Tumor');
        ylabel('[VEGFR2]_t (pM)');
        xlabel('Time (days)');

        subplot(5,3,12);
        hold on;
        plot(t_days, VEGFR2_r_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFR2_r_all,2),'k','LineWidth',lw_mean);
        title('VEGFR2 – Rest of Body');
        ylabel('[VEGFR2]_r (pM)');
        xlabel('Time (days)');

        % ---- Row 5: VEGF–VEGFR2 (pM) ----
        subplot(5,3,13);
        hold on;
        plot(t_days, VEGFVEGFR2_b_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFVEGFR2_b_all,2),'k','LineWidth',lw_mean);
        title('VEGF–VEGFR2 – Blood');
        ylabel('[VEGF–VEGFR2]_b (pM)');
        xlabel('Time (days)');

        subplot(5,3,14);
        hold on;
        plot(t_days, VEGFVEGFR2_t_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFVEGFR2_t_all,2),'k','LineWidth',lw_mean);
        title('VEGF–VEGFR2 – Tumor');
        ylabel('[VEGF–VEGFR2]_t (pM)');
        xlabel('Time (days)');

        subplot(5,3,15);
        hold on;
        plot(t_days, VEGFVEGFR2_r_all,'LineWidth',lw_ind);
        plot(t_days, mean(VEGFVEGFR2_r_all,2),'k','LineWidth',lw_mean);
        title('VEGF–VEGFR2 – Rest of Body');
        ylabel('[VEGF–VEGFR2]_r (pM)');
        xlabel('Time (days)');

        sgtitle(['Task 3 - PopPK ', scenarios(s).name]);

        saveas(fig, scenarios(s).filename);
        close(fig);
    end
end

%% Helper: simulate given regimen for one patient
function [T_all, Y_all] = simulate_regimen(regimen, p_local, y_ss, ...
                                           solver, options, ...
                                           dose10_conc_pM, dose1_conc_pM, ...
                                           tspan_single, ...
                                           metro_dose_days, metro_total_days, ...
                                           minutes_per_day, n)
    switch lower(regimen)
        case 'single'
            % Single 10 mg/kg bolus at t = 0, simulate 21 days
            y0 = y_ss;
            y0(n.Ab_b) = y0(n.Ab_b) + dose10_conc_pM;
            [T_all, Y_all] = solver(@VEGFAbeqns, tspan_single, y0, options, p_local, n);

        case 'metronomic'
            % 1 mg/kg equivalent per day x metro_dose_days, total metro_total_days
            y_curr = y_ss;
            t_curr = 0;

            T_all = [];
            Y_all = [];

            % Dosing phase
            for d = 1:metro_dose_days
                % dose for this day
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

            % Tail phase after last dose
            final_time = metro_total_days * minutes_per_day;
            if t_curr < final_time
                tspan_tail = t_curr:1:final_time;
                [T_tail, Y_tail] = solver(@VEGFAbeqns, tspan_tail, y_curr, options, p_local, n);
                T_all = [T_all; T_tail(2:end)];
                Y_all = [Y_all; Y_tail(2:end,:)];
            end

        otherwise
            error('Unknown regimen: %s', regimen);
    end
end

%% Helper: sample from N(mu,sigma) but truncated at > 0
function vals = sample_positive_normal(mu,sigma,N)
    vals = zeros(N,1);
    for i = 1:N
        v = -1;
        while v <= 0
            v = normrnd(mu, sigma);
        end
        vals(i) = v;
    end
end
