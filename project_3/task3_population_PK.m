function task3_population_PK

    clearvars -except task3_population_PK;
    clc;

    %======================================================================
    % Task 3 – Population PK: virtual patients with variable clearance
    % and rest-of-body volume
    %
    % Based on the simplified Stefanini et al. model:
    % one VEGF isoform, one receptor, one antibody
    % three tissues (blood, tumor, rest of body)
    %
    % This script creates a virtual population of patients by varying:
    %   - kcl_A (Ab clearance) and kcl_VA (VEGF–Ab clearance), CV = 50%
    %   - Vol_r (rest-of-body volume), CV = 25%
    % and simulates a single 10 mg/kg IV bolus of bevacizumab.
    %
    % It then visualizes VEGF, Ab, and VEGF–Ab in:
    %   - blood
    %   - tumor
    %   - rest-of-body ("normal tissue")
    %
    % Professor suggestions implemented:
    %   - Use ode15s instead of ode45
    %   - Increase MaxStep to 5e-1
    %   - Use parfor to parallelize patient loop
    %======================================================================

    %-----------------------------
    % Index array (for readability)
    %-----------------------------
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

    %-----------------------------------------------
    % Parameters – compartment volumes (interstitial)
    %-----------------------------------------------
    p.Vol_b = 5 * 0.60;    % blood: 5 L, 60% available
    p.Vol_t = 1 * 0.611;   % tumor: 1 L, 61.1% available
    p.Vol_r = 40 * 0.0816; % rest of body: 40 L, 8.16% available

    %-------------------------------
    % Parameters – clearance, intnl
    %-------------------------------
    p.kcl_V   = 0.0648;   % VEGF clearance (min^-1)
    p.kcl_VA  = 2.2e-5;   % VEGF–Ab complex clearance (min^-1)
    p.kcl_A   = 2.2e-5;   % Ab clearance (min^-1)

    p.kintR2  = 0.0168;   % receptor internalization (min^-1)
    p.kintVR2 = 0.0168;   % receptor–VEGF internalization (min^-1)

    %-----------------------
    % Parameters – production
    %-----------------------
    p.VEGFprod_b   = 0;
    p.VEGFprod_t   = 3;
    p.VEGFprod_r   = 10;

    p.VEGFR2prod_b = 0;
    p.VEGFR2prod_t = 2.85e2 * p.kintR2;
    p.VEGFR2prod_r = 2.2e3  * p.kintR2;

    %------------------------
    % Parameters – binding
    %------------------------
    p.konVR  = 6e-4;    % VEGF–VEGFR2 on-rate (pM^-1 min^-1)
    p.koffVR = 0.06;    % VEGF–VEGFR2 off-rate (min^-1)

    p.konVA  = 5.52e-6; % VEGF–Ab on-rate (pM^-1 min^-1)
    p.koffVA = 0.012;   % VEGF–Ab off-rate (min^-1)

    %-----------------------------
    % Compartmental transport
    %-----------------------------
    p.k_bt = 4.12e-4;           % blood -> tumor
    p.k_br = 3.19e-3;           % blood -> rest of body
    p.k_tb = 4.12e-4;           % tumor -> blood
    p.k_rb = 3.19e-3 + 6.15e-4; % rest of body -> blood (includes lymphatics)

    % AbEx flag: 1 = normal extravasation, 0 = no transport
    p.AbEx = 1;

    %-----------------------------
    % Antibody dose: 10 mg/kg IV
    %-----------------------------
    dose_mg_total   = 10 * 70;         % 10 mg/kg * 70 kg
    MW_mg_per_mol   = 150000000;       % 150 kDa in mg/mol
    dose_mol_total  = dose_mg_total / MW_mg_per_mol;
    dose_pmol_total = dose_mol_total * 1e12; % pmol
    dose_conc_pM    = dose_pmol_total / p.Vol_b; % pM in blood

    %-----------------------------
    % ODE solver (stiff) + options
    %-----------------------------
    solver = @ode15s;  % stiff solver

    options = odeset('MaxStep', 5e-1, ...   % larger MaxStep (prof suggestion)
                     'AbsTol',  1e-5, ...
                     'RelTol',  1e-5, ...
                     'InitialStep', 1e-2);

    %======================================================================
    % Population definition
    %======================================================================
    N_patients = 100;

    % Baseline (mean) values
    base_kcl_A  = p.kcl_A;
    base_kcl_VA = p.kcl_VA;
    base_Vol_r  = p.Vol_r;

    CV_kcl = 0.50;   % 50% coeff. of variation for clearances
    CV_Vr  = 0.25;   % 25% coeff. of variation for rest-of-body volume

    % For reproducibility
    rng(1);

    % Sample individual parameter values (truncate at > 0, via helper)
    kcl_A_samples  = sample_positive_normal(base_kcl_A,  CV_kcl * base_kcl_A,  N_patients);
    kcl_VA_samples = sample_positive_normal(base_kcl_VA, CV_kcl * base_kcl_VA, N_patients);
    Vol_r_samples  = sample_positive_normal(base_Vol_r,  CV_Vr  * base_Vol_r,  N_patients);

    %======================================================================
    % Simulation settings
    %======================================================================
    % Steady-state (no Ab) run time
    sstime  = 60 * 24 * 10;     % 10 days, minutes
    tspanSS = 0:1:sstime;

    % Post-dose simulation time
    endtime = 60 * 24 * 21;     % 21 days post-dose
    tspan   = 0:1:endtime;
    nt      = numel(tspan);

    % Preallocate outputs for blood, tumor, rest-of-body
    VEGF_b_all     = zeros(nt, N_patients);
    Ab_b_all       = zeros(nt, N_patients);
    VEGFAb_b_all   = zeros(nt, N_patients);

    VEGF_t_all     = zeros(nt, N_patients);
    Ab_t_all       = zeros(nt, N_patients);
    VEGFAb_t_all   = zeros(nt, N_patients);

    VEGF_r_all     = zeros(nt, N_patients);
    Ab_r_all       = zeros(nt, N_patients);
    VEGFAb_r_all   = zeros(nt, N_patients);

    % Time in days for plotting
    t_days = tspan / (60 * 24);

    %======================================================================
    % Loop over virtual patients (PARALLEL)
    %======================================================================
    % If parfor is problematic on your machine, you can safely replace
    % "parfor" with "for".
    parfor ip = 1:N_patients

        %-------------------------------------
        % Patient-specific parameter struct
        %-------------------------------------
        p_i = p;
        p_i.kcl_A  = kcl_A_samples(ip);
        p_i.kcl_VA = kcl_VA_samples(ip);
        p_i.Vol_r  = Vol_r_samples(ip);

        %-------------------------------------
        % 1) Run to steady state with no Ab
        %-------------------------------------
        y0_SS = zeros(15,1);
        [~, Yss] = solver(@VEGFAbeqns, tspanSS, y0_SS, options, p_i, n);
        y_ss = Yss(end, :).';   % steady-state concentrations

        %-------------------------------------
        % 2) Add single IV bolus of Ab to blood
        %-------------------------------------
        y0_post = y_ss;
        y0_post(n.Ab_b) = y0_post(n.Ab_b) + dose_conc_pM;

        %-------------------------------------
        % 3) Simulate 21 days post-dose
        %-------------------------------------
        [~, Y_i] = solver(@VEGFAbeqns, tspan, y0_post, options, p_i, n);

        % Store concentrations in all 3 compartments
        VEGF_b_all(:, ip)   = Y_i(:, n.VEGF_b);
        Ab_b_all(:, ip)     = Y_i(:, n.Ab_b);
        VEGFAb_b_all(:, ip) = Y_i(:, n.VEGFAb_b);

        VEGF_t_all(:, ip)   = Y_i(:, n.VEGF_t);
        Ab_t_all(:, ip)     = Y_i(:, n.Ab_t);
        VEGFAb_t_all(:, ip) = Y_i(:, n.VEGFAb_t);

        VEGF_r_all(:, ip)   = Y_i(:, n.VEGF_r);
        Ab_r_all(:, ip)     = Y_i(:, n.Ab_r);
        VEGFAb_r_all(:, ip) = Y_i(:, n.VEGFAb_r);
    end

    %======================================================================
    % Visualization: 3 × 3 figure
    % Rows: VEGF, Ab, VEGF–Ab
    % Cols: Blood, Tumor, Rest-of-body
    %======================================================================
    fig = figure('Visible','off');
    set(fig, 'Position', [100 100 1400 900]);

    lw_ind  = 0.4;   % individual curves
    lw_mean = 2;     % population mean

    %-----------------------------
    % Row 1: VEGF
    %-----------------------------
    % Blood
    subplot(3,3,1);
    hold on;
    plot(t_days, VEGF_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_b_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Blood');
    xlabel('Time (days)');
    ylabel('[VEGF]_b (pM)');

    % Tumor
    subplot(3,3,2);
    hold on;
    plot(t_days, VEGF_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_t_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Tumor');
    xlabel('Time (days)');
    ylabel('[VEGF]_t (pM)');

    % Rest-of-body
    subplot(3,3,3);
    hold on;
    plot(t_days, VEGF_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGF_r_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF – Rest of body');
    xlabel('Time (days)');
    ylabel('[VEGF]_r (pM)');

    %-----------------------------
    % Row 2: Free Ab
    %-----------------------------
    % Blood
    subplot(3,3,4);
    hold on;
    plot(t_days, Ab_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(Ab_b_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('Ab – Blood (PK)');
    xlabel('Time (days)');
    ylabel('[Ab]_b (pM)');

    % Tumor
    subplot(3,3,5);
    hold on;
    plot(t_days, Ab_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(Ab_t_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('Ab – Tumor');
    xlabel('Time (days)');
    ylabel('[Ab]_t (pM)');

    % Rest-of-body
    subplot(3,3,6);
    hold on;
    plot(t_days, Ab_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(Ab_r_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('Ab – Rest of body');
    xlabel('Time (days)');
    ylabel('[Ab]_r (pM)');

    %-----------------------------
    % Row 3: VEGF–Ab complex
    %-----------------------------
    % Blood
    subplot(3,3,7);
    hold on;
    plot(t_days, VEGFAb_b_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFAb_b_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–Ab – Blood');
    xlabel('Time (days)');
    ylabel('[VEGF–Ab]_b (pM)');

    % Tumor
    subplot(3,3,8);
    hold on;
    plot(t_days, VEGFAb_t_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFAb_t_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–Ab – Tumor');
    xlabel('Time (days)');
    ylabel('[VEGF–Ab]_t (pM)');

    % Rest-of-body
    subplot(3,3,9);
    hold on;
    plot(t_days, VEGFAb_r_all, 'LineWidth', lw_ind);
    plot(t_days, mean(VEGFAb_r_all, 2), 'k', 'LineWidth', lw_mean);
    hold off;
    title('VEGF–Ab – Rest of body');
    xlabel('Time (days)');
    ylabel('[VEGF–Ab]_r (pM)');

    sgtitle('Task 3 – Population PK: VEGF, Ab, and VEGF–Ab in Blood, Tumor, and Normal Tissue');

    saveas(fig, 'task3_population_PK.png');
    close(fig);

end

%======================================================================
% Helper function: sample from N(mu, sigma) but ensure strictly positive
%======================================================================
function vals = sample_positive_normal(mu, sigma, N)
    vals = zeros(N,1);
    for i = 1:N
        v = -1;
        while v <= 0
            v = normrnd(mu, sigma);
        end
        vals(i) = v;
    end
end
