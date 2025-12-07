function task1_single_injection

clearvars -except task1_single_injection;
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

dose_mg_total   = 10 * 70;
MW_mg_per_mol   = 150000000;
dose_mol_total  = dose_mg_total / MW_mg_per_mol;
dose_pmol_total = dose_mol_total * 1e12;
dose_conc_pM    = dose_pmol_total / p.Vol_b;

options = odeset('MaxStep', 5e-2, ...
                 'AbsTol', 1e-5, ...
                 'RelTol', 1e-5, ...
                 'InitialStep', 1e-2);

y0      = zeros(15,1);
sstime  = 60 * 24 * 10;
tspanSS = 0:1:sstime;

[~, Yss] = ode45(@VEGFAbeqns, tspanSS, y0, options, p, n);
y_ss = Yss(end, :).';

endtime = 60 * 24 * 21;
tspan   = 0:1:endtime;

p_noEx          = p;
p_noEx.AbEx     = 0;
y0_noEx         = y_ss;
y0_noEx(n.Ab_b) = y0_noEx(n.Ab_b) + dose_conc_pM;
[T_noEx, Y_noEx] = ode45(@VEGFAbeqns, tspan, y0_noEx, options, p_noEx, n);

p_Ex          = p;
p_Ex.AbEx     = 1;
y0_Ex         = y_ss;
y0_Ex(n.Ab_b) = y0_Ex(n.Ab_b) + dose_conc_pM;
[T_Ex, Y_Ex] = ode45(@VEGFAbeqns, tspan, y0_Ex, options, p_Ex, n);

t_noEx_days = T_noEx / (60*24);
t_Ex_days   = T_Ex   / (60*24);

Ab_noEx_uM_r = Y_noEx(:,n.Ab_r) / 1e6;
Ab_noEx_uM_b = Y_noEx(:,n.Ab_b) / 1e6;
Ab_noEx_uM_t = Y_noEx(:,n.Ab_t) / 1e6;

Ab_Ex_uM_r = Y_Ex(:,n.Ab_r) / 1e6;
Ab_Ex_uM_b = Y_Ex(:,n.Ab_b) / 1e6;
Ab_Ex_uM_t = Y_Ex(:,n.Ab_t) / 1e6;

C_noEx_nM_r = Y_noEx(:,n.VEGFAb_r) / 1e3;
C_noEx_nM_b = Y_noEx(:,n.VEGFAb_b) / 1e3;
C_noEx_nM_t = Y_noEx(:,n.VEGFAb_t) / 1e3;

C_Ex_nM_r = Y_Ex(:,n.VEGFAb_r) / 1e3;
C_Ex_nM_b = Y_Ex(:,n.VEGFAb_b) / 1e3;
C_Ex_nM_t = Y_Ex(:,n.VEGFAb_t) / 1e3;

ymax_complex_Ex = max(C_Ex_nM_t) * 1.1;

fig = figure('Visible','off');
set(fig, 'Position', [100 100 1200 800]);
lw = 2;

% Row 1: Free VEGF (pM)
% --- FIX: zoom y-limits so blood 0–5 pM drop/recovery is visible ---
subplot(3,2,1);
plot(t_noEx_days, Y_noEx(:,n.VEGF_r), '-',  ...
     t_noEx_days, Y_noEx(:,n.VEGF_b), '--', ...
     t_noEx_days, Y_noEx(:,n.VEGF_t), ':',  'LineWidth', lw);
title('Free VEGF – NO Ab extravasation');
xlabel('Time (days)'); ylabel('Free VEGF (pM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');
ylim([0 220]);

subplot(3,2,2);
plot(t_Ex_days, Y_Ex(:,n.VEGF_r), '-', ...
     t_Ex_days, Y_Ex(:,n.VEGF_b), '--', ...
     t_Ex_days, Y_Ex(:,n.VEGF_t), ':', 'LineWidth', lw);
title('Free VEGF – WITH Ab extravasation');
xlabel('Time (days)'); ylabel('[VEGF] (pM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');
ylim([0 220]);

subplot(3,2,3);
plot(t_noEx_days, Ab_noEx_uM_r, '-', ...
     t_noEx_days, Ab_noEx_uM_b, '--', ...
     t_noEx_days, Ab_noEx_uM_t, ':', 'LineWidth', lw);
title('Free anti-VEGF – NO Ab extravasation');
xlabel('Time (days)'); ylabel('Free anti-VEGF (\muM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');

subplot(3,2,4);
plot(t_Ex_days, Ab_Ex_uM_r, '-', ...
     t_Ex_days, Ab_Ex_uM_b, '--', ...
     t_Ex_days, Ab_Ex_uM_t, ':', 'LineWidth', lw);
title('Free anti-VEGF – WITH Ab extravasation');
xlabel('Time (days)'); ylabel('Free anti-VEGF (\muM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');

subplot(3,2,5);
plot(t_noEx_days, C_noEx_nM_r, '-', ...
     t_noEx_days, C_noEx_nM_b, '--', ...
     t_noEx_days, C_noEx_nM_t, ':', 'LineWidth', lw);
title('VEGF–anti-VEGF complex – NO Ab extravasation');
xlabel('Time (days)'); ylabel('[VEGF–anti-VEGF] (nM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');

subplot(3,2,6);
plot(t_Ex_days, C_Ex_nM_r, '-', ...
     t_Ex_days, C_Ex_nM_b, '--', ...
     t_Ex_days, C_Ex_nM_t, ':', 'LineWidth', lw);
title('VEGF–anti-VEGF complex – WITH Ab extravasation');
xlabel('Time (days)'); ylabel('[VEGF–anti-VEGF] (nM)');
legend({'Normal tissue','Blood','Tumor'}, 'Location','best');
ylim([0 ymax_complex_Ex]);

sgtitle('Task 1 – Single 10 mg/kg IV injection of bevacizumab');
saveas(fig, 'task1_single_injection.png');
close(fig);
end
