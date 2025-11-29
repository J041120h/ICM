function task1_single_injection
% Task 1 – Single 10 mg/kg IV injection of bevacizumab
% Uses VEGFAbeqns.m (must be in same folder, UNMODIFIED).
% Produces a 3x2 panel figure and saves it as task1_single_injection.png.
%
% Panels:
%   Row 1: Free VEGF (no extravasation vs with extravasation)
%   Row 2: Free antibody
%   Row 3: VEGF–antibody complex

    clearvars -except task1_single_injection;  %#ok<CLVAR>
    clc;

    %% Index array (must match VEGFAbeqns.m)
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

    %% Parameters (from starter code / assignment)
    % Compartment volumes (L)
    p.Vol_b = 5 * 0.60;    % blood plasma
    p.Vol_t = 1 * 0.611;   % tumor interstitial
    p.Vol_r = 40 * 0.0816; % rest-of-body interstitial

    % Clearance and internalization (min^-1)
    p.kcl_V   = 0.0648;    % free VEGF clearance
    p.kcl_VA  = 2.2e-5;    % VEGF–Ab clearance
    p.kcl_A   = 2.2e-5;    % Ab clearance

    p.kintR2  = 0.0168;    % VEGFR2 internalization
    p.kintVR2 = 0.0168;    % VEGF–VEGFR2 internalization

    % Production (pM/min)
    p.VEGFprod_b   = 0;
    p.VEGFprod_t   = 3;
    p.VEGFprod_r   = 10;

    p.VEGFR2prod_b = 0;
    p.VEGFR2prod_t = 2.85e2 * p.kintR2;   % 285 pM * kint
    p.VEGFR2prod_r = 2.2e3  * p.kintR2;   % 2200 pM * kint

    % Binding kinetics (kon: pM^-1 min^-1, koff: min^-1)
    p.konVR  = 6e-4;
    p.koffVR = 0.06;

    p.konVA  = 5.52e-6;
    p.koffVA = 0.012;

    % Transport rates (min^-1)
    p.k_bt = 4.12e-4;               % blood -> tumor
    p.k_br = 3.19e-3;               % blood -> rest
    p.k_tb = 4.12e-4;               % tumor -> blood
    p.k_rb = 3.19e-3 + 6.15e-4;     % rest -> blood (includes lymphatics)

    % Default: Ab transport allowed (will override for each case)
    p.AbEx = 1;

    %% Antibody dose (single IV bolus)
    % 10 mg/kg * 70 kg = 700 mg total; MW = 150 kDa = 150,000,000 mg/mol
    dose_mg_total   = 10 * 70;
    MW_mg_per_mol   = 150000000;
    dose_mol_total  = dose_mg_total / MW_mg_per_mol;
    dose_pmol_total = dose_mol_total * 1e12;    % pmol
    dose_conc_pM    = dose_pmol_total / p.Vol_b;  % pM increase in blood

    %% ODE solver options
    options = odeset('MaxStep', 5e-2, ...
                     'AbsTol', 1e-5, ...
                     'RelTol', 1e-5, ...
                     'InitialStep', 1e-2);

    %% 1) Run to steady state with NO antibody
    y0      = zeros(15,1);
    sstime  = 60 * 24 * 10;          % 10 days in minutes
    tspanSS = 0:1:sstime;

    % Antibody-related parameters irrelevant here (Ab=0)
    [Tss, Yss] = ode45(@VEGFAbeqns, tspanSS, y0, options, p, n); %#ok<ASGLU>
    y_ss = Yss(end, :).';   % steady-state vector

    %% 2) Post-dose simulations (21 days after dose)
    endtime = 60 * 24 * 21;          % 21 days in minutes
    tspan   = 0:1:endtime;

    % ---------- Case A: NO antibody extravasation ----------
    p_noEx       = p;
    p_noEx.AbEx  = 0;                       % key: no extravasation
    y0_noEx      = y_ss;
    y0_noEx(n.Ab_b) = y0_noEx(n.Ab_b) + dose_conc_pM;

    [T_noEx, Y_noEx] = ode45(@VEGFAbeqns, tspan, y0_noEx, options, p_noEx, n);

    % ---------- Case B: WITH antibody extravasation ----------
    p_Ex       = p;
    p_Ex.AbEx  = 1;                         % key: extravasation allowed
    y0_Ex      = y_ss;
    y0_Ex(n.Ab_b) = y0_Ex(n.Ab_b) + dose_conc_pM;

    [T_Ex, Y_Ex] = ode45(@VEGFAbeqns, tspan, y0_Ex, options, p_Ex, n);

    % Convert time to days for plotting
    t_noEx_days = T_noEx / (60*24);
    t_Ex_days   = T_Ex   / (60*24);

    %% 3) Plot (but do NOT show window) and save locally
    fig = figure('Visible','off');
    set(fig, 'Position', [100 100 1200 800]);

    lw = 2;

    % Row 1: Free VEGF
    subplot(3,2,1);
    plot(t_noEx_days, Y_noEx(:,n.VEGF_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.VEGF_b), '--', ...
         t_noEx_days, Y_noEx(:,n.VEGF_t), ':',  'LineWidth', lw);
    title('Free VEGF – NO Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,2);
    plot(t_Ex_days, Y_Ex(:,n.VEGF_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.VEGF_b), '--', ...
         t_Ex_days, Y_Ex(:,n.VEGF_t), ':',  'LineWidth', lw);
    title('Free VEGF – WITH Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    % Row 2: Free antibody
    subplot(3,2,3);
    plot(t_noEx_days, Y_noEx(:,n.Ab_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.Ab_b), '--', ...
         t_noEx_days, Y_noEx(:,n.Ab_t), ':',  'LineWidth', lw);
    title('Free antibody – NO Ab extravasation');
    xlabel('Time (days)'); ylabel('[Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,4);
    plot(t_Ex_days, Y_Ex(:,n.Ab_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.Ab_b), '--', ...
         t_Ex_days, Y_Ex(:,n.Ab_t), ':',  'LineWidth', lw);
    title('Free antibody – WITH Ab extravasation');
    xlabel('Time (days)'); ylabel('[Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    % Row 3: VEGF–Ab complex
    subplot(3,2,5);
    plot(t_noEx_days, Y_noEx(:,n.VEGFAb_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.VEGFAb_b), '--', ...
         t_noEx_days, Y_noEx(:,n.VEGFAb_t), ':',  'LineWidth', lw);
    title('VEGF–Ab complex – NO Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF–Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,6);
    plot(t_Ex_days, Y_Ex(:,n.VEGFAb_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.VEGFAb_b), '--', ...
         t_Ex_days, Y_Ex(:,n.VEGFAb_t), ':',  'LineWidth', lw);
    title('VEGF–Ab complex – WITH Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF–Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    sgtitle('Task 1 – Single 10 mg/kg IV injection of bevacizumab');

    % Save figure locally and close
    saveas(fig, 'task1_single_injection.png');
    close(fig);
end