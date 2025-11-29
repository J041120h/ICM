function task2_metronomic

    clearvars -except task2_metronomic;
    clc;

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

    p.AbEx = 1;

    dose10_mg_total   = 10 * 70;
    MW_mg_per_mol     = 150000000;
    dose10_mol_total  = dose10_mg_total / MW_mg_per_mol;
    dose10_pmol_total = dose10_mol_total * 1e12;
    dose10_conc_pM    = dose10_pmol_total / p.Vol_b;

    dose1_conc_pM = dose10_conc_pM / 10;

    options = odeset('MaxStep', 5e-2, ...
                     'AbsTol', 1e-5, ...
                     'RelTol', 1e-5, ...
                     'InitialStep', 1e-2);

    y0      = zeros(15,1);
    sstime  = 60 * 24 * 10;
    tspanSS = 0:1:sstime;

    [~, Yss] = ode45(@VEGFAbeqns, tspanSS, y0, options, p, n);
    y_ss = Yss(end, :).';

    n_days_dose = 10;
    total_days  = 30;
    minutes_per_day = 24 * 60;

    function [T_all, Y_all] = run_metronomic_for_p(p_local)
        y_curr  = y_ss;
        t_curr  = 0;
        T_all   = [];
        Y_all   = [];
        for d = 1:n_days_dose
            y_curr(n.Ab_b) = y_curr(n.Ab_b) + dose1_conc_pM;
            tspan_day = t_curr:1:(t_curr + minutes_per_day);
            [T_day, Y_day] = ode45(@VEGFAbeqns, tspan_day, y_curr, options, p_local, n);
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
        final_time = total_days * minutes_per_day;
        if t_curr < final_time
            tspan_tail = t_curr:1:final_time;
            [T_tail, Y_tail] = ode45(@VEGFAbeqns, tspan_tail, y_curr, options, p_local, n);
            T_all = [T_all; T_tail(2:end)];
            Y_all = [Y_all; Y_tail(2:end,:)];
        end
    end

    p_noEx      = p;
    p_noEx.AbEx = 0;

    [T_noEx, Y_noEx] = run_metronomic_for_p(p_noEx);

    p_Ex      = p;
    p_Ex.AbEx = 1;

    [T_Ex, Y_Ex] = run_metronomic_for_p(p_Ex);

    t_noEx_days = T_noEx / (60*24);
    t_Ex_days   = T_Ex   / (60*24);

    fig = figure('Visible','off');
    set(fig, 'Position', [100 100 1200 800]);
    lw = 2;

    subplot(3,2,1);
    plot(t_noEx_days, Y_noEx(:,n.VEGF_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.VEGF_b), '--', ...
         t_noEx_days, Y_noEx(:,n.VEGF_t), ':',  'LineWidth', lw);
    title('Free VEGF – Metronomic, NO Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,2);
    plot(t_Ex_days, Y_Ex(:,n.VEGF_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.VEGF_b), '--', ...
         t_Ex_days, Y_Ex(:,n.VEGF_t), ':',  'LineWidth', lw);
    title('Free VEGF – Metronomic, WITH Ab extravasation');
    xlabel('Time (days)'); ylabel('[VEGF] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,3);
    plot(t_noEx_days, Y_noEx(:,n.Ab_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.Ab_b), '--', ...
         t_noEx_days, Y_noEx(:,n.Ab_t), ':',  'LineWidth', lw);
    title('Free Ab – Metronomic, NO extravasation');
    xlabel('Time (days)'); ylabel('[Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,4);
    plot(t_Ex_days, Y_Ex(:,n.Ab_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.Ab_b), '--', ...
         t_Ex_days, Y_Ex(:,n.Ab_t), ':',  'LineWidth', lw);
    title('Free Ab – Metronomic, WITH extravasation');
    xlabel('Time (days)'); ylabel('[Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,5);
    plot(t_noEx_days, Y_noEx(:,n.VEGFAb_r), '-',  ...
         t_noEx_days, Y_noEx(:,n.VEGFAb_b), '--', ...
         t_noEx_days, Y_noEx(:,n.VEGFAb_t), ':',  'LineWidth', lw);
    title('VEGF–Ab – Metronomic, NO extravasation');
    xlabel('Time (days)'); ylabel('[VEGF–Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    subplot(3,2,6);
    plot(t_Ex_days, Y_Ex(:,n.VEGFAb_r), '-',  ...
         t_Ex_days, Y_Ex(:,n.VEGFAb_b), '--', ...
         t_Ex_days, Y_Ex(:,n.VEGFAb_t), ':',  'LineWidth', lw);
    title('VEGF–Ab – Metronomic, WITH extravasation');
    xlabel('Time (days)'); ylabel('[VEGF–Ab] (pM)');
    legend({'Rest of body','Blood','Tumor'}, 'Location','best');

    sgtitle('Task 2 – Metronomic therapy (1 mg/kg daily × 10 days)');

    saveas(fig, 'task2_metronomic.png');
    close(fig);
end
