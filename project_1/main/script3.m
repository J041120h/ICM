% ========================================================================
% Project 1 – Estimating Cardiac Output from ABP 
% Team 1: Daniel Chen, Harry Jiang, Carol Li, Amanda Xu
% Problems 1, 3, 5 (with C1/C2 calibration fix) and 6
% ========================================================================

%% -------------------------- PROBLEM 1 ----------------------------------
clear; clc;
fprintf('\n=== Starting Problem 1: Producing ABP traces for Patient #20 ===\n');

PATH_ABP = 's00020-2567-03-30-17-47_ABP.txt';
PATH_NUM = 's00020-2567-03-30-17-47n.txt';
OUT_DIR  = '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/result';

if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

T10 = 10*3600; T11 = 11*3600;

raw = readmatrix(PATH_ABP,'FileType','text');
t   = raw(:,1); abp = raw(:,2);
t = t(:); abp = abp(:);

OnsetTimes = wabp(abp);

start10_idx = find(t >= T10, 1, 'first');
k10 = find(OnsetTimes >= start10_idx, 1, 'first');
onsets10 = OnsetTimes(k10:(k10+20));

start11_idx = find(t >= T11, 1, 'first');
k11 = find(OnsetTimes >= start11_idx, 1, 'first');
onsets11 = OnsetTimes(k11:(k11+20));

feat10 = abpfeature(abp, onsets10); 
feat11 = abpfeature(abp, onsets11);

[BeatQ10, ~] = jSQI(feat10, onsets10(1:end-1), abp);
[BeatQ11, ~] = jSQI(feat11, onsets11(1:end-1), abp);

fprintf('\nPlotting traces for 10h and 11h\n');

% --- 10h plot ---
i0 = onsets10(1); i1 = onsets10(end);
tt = t(i0:i1); yy = abp(i0:i1);
on20 = onsets10(1:end-1); eos03 = feat10(:,9); eosMin = feat10(:,11);

f1 = figure('Color','w');
plot(tt,yy,'LineWidth',1); hold on;
plot(t(on20),abp(on20),'*','MarkerSize',5);
plot(t(eos03),abp(eos03),'xr','MarkerSize',7,'LineWidth',1.2);
plot(t(eosMin),abp(eosMin),'o','MarkerSize',5);
xlabel('Time (s)'); ylabel('ABP (mmHg)');
title('Patient #20, Starting at 10 h');
box on; grid on;
print(fullfile(OUT_DIR,'Problem1_10h.png'),'-dpng','-r300'); close(gcf);

% --- 11h plot ---
i0 = onsets11(1); i1 = onsets11(end);
tt = t(i0:i1); yy = abp(i0:i1);
on20 = onsets11(1:end-1); eos03 = feat11(:,9); eosMin = feat11(:,11);

f2 = figure('Color','w');
plot(tt,yy,'LineWidth',1); hold on;
plot(t(on20),abp(on20),'*','MarkerSize',5);
plot(t(eos03),abp(eos03),'xr','MarkerSize',7,'LineWidth',1.2);
plot(t(eosMin),abp(eosMin),'o','MarkerSize',5);
xlabel('Time (s)'); ylabel('ABP (mmHg)');
title('Patient #20, Starting at 11 h');
box on; grid on;
ax=findall(gcf,'type','axes'); set(ax,'Toolbar',[]);
print(fullfile(OUT_DIR,'Problem1_11h.png'),'-dpng','-r300'); close(gcf);
fprintf('Saved Figure to %s\n',OUT_DIR);
fprintf('\n=== Problem 1 Completed ===\n');


%% -------------------------- PROBLEM 3 ----------------------------------
fprintf('\n=== Starting Problem 3: Continuous CO estimation (Liljestrand) ===\n');
Fs=125; TMAX_HR=12; TMAX_SEC=TMAX_HR*3600;
mask12=(t<=TMAX_SEC); time12=t(mask12); ABP12=abp(mask12);

t_on=wabp(ABP12);
feat=abpfeature(ABP12,t_on);
beatq=jSQI(feat,t_on(1:end-1),ABP12);
time=time12; ABP=ABP12;
save('patient20_first12h.mat','time','ABP','t_on','feat','beatq');

[co_uncal,to_min,~,fea]=estimateCO_v2('patient20_first12h.mat',5,1);
to_hr=to_min/60;
Psys=fea(:,2); Pdias=fea(:,4); PP=fea(:,5); MAP=fea(:,6);
Period=fea(:,7); HR=60*Fs./Period;

tbl=readtable(PATH_NUM,'FileType','text','HeaderLines',2);
time_num_all=tbl{:,1}; CO_TD_all=tbl{:,end};
nz_idx_all=find(CO_TD_all~=0 & ~isnan(CO_TD_all));
td_time_sec_all=time_num_all(nz_idx_all);
td_CO_Lmin_all=CO_TD_all(nz_idx_all);
mask12_TD=td_time_sec_all<=TMAX_SEC;
td_time_sec_12h=td_time_sec_all(mask12_TD);
td_CO_Lmin_12h=td_CO_Lmin_all(mask12_TD);
td_time_hr_12h=td_time_sec_12h/3600;

TD0=td_CO_Lmin_all(1); t0_sec=td_time_sec_all(1);
[~,i_near]=min(abs(to_min-(t0_sec/60)));
k_cal=TD0/co_uncal(i_near);
co_cal=k_cal*co_uncal;

PP_td=[]; MAP_td=[]; HR_td=[];
for ii=1:numel(td_time_sec_12h)
    [~,jj]=min(abs(to_min-(td_time_sec_12h(ii)/60)));
    PP_td(ii,1)=PP(jj); MAP_td(ii,1)=MAP(jj); HR_td(ii,1)=HR(jj);
end

figure('Color','w','Name','Problem3 Liljestrand');
sgtitle('Problem 3 – Patient #20 | Liljestrand (#5) | First 12h');
subplot(4,1,1); hold on;
plot(to_hr,co_cal,'Color',[0.5 0.5 0.5],'LineWidth',1.3);
stem(td_time_hr_12h,td_CO_Lmin_12h,'k','filled');
ylabel('CO [L/min]'); xlim([0 TMAX_HR]); grid on;
title(sprintf('Calibrated CO (k = %.3f)',k_cal));
subplot(4,1,2); hold on;
plot(to_hr,PP,'Color',[0.5 0.5 0.5]); stem(td_time_hr_12h,PP_td,'k','filled');
ylabel('PP [mmHg]'); xlim([0 TMAX_HR]); grid on;
subplot(4,1,3); hold on;
plot(to_hr,MAP,'Color',[0.5 0.5 0.5]); stem(td_time_hr_12h,MAP_td,'k','filled');
ylabel('MAP [mmHg]'); xlim([0 TMAX_HR]); grid on;
subplot(4,1,4); hold on;
plot(to_hr,HR,'Color',[0.5 0.5 0.5]); stem(td_time_hr_12h,HR_td,'k','filled');
ylabel('HR [bpm]'); xlabel('Time [hours]');
xlim([0 TMAX_HR]); grid on;
ax=findall(gcf,'type','axes'); set(ax,'Toolbar',[]);
saveas(gcf,fullfile(OUT_DIR,'Problem3_patient20.png'));
fprintf('\n=== Problem 3 Completed ===\n');


%% -------------------------- PROBLEM 5 ----------------------------------
fprintf('\n=== Starting Problem 5: Parlikar (#14) with C1 calibration ===\n');
ABP_FILES = {'s00020-2567-03-30-17-47_ABP.txt','s10013-2564-11-01-23-37_ABP.txt','s10205-2631-06-14-11-38_ABP.txt','s10241-3150-11-15-13-27_ABP.txt'};
NUM_FILES = {'s00020-2567-03-30-17-47n.txt','s10013-2564-11-01-23-37n.txt','s10205-2631-06-14-11-38n.txt','s10241-3150-11-15-13-27n.txt'};
OUT_ROOT  = '/Users/harry/Desktop/ICM/Project 1/PhysioToolkitCardiacOutput_MatlabCode/result';

Fs=125; TMAX_HR=12; TMAX_SEC=TMAX_HR*3600;
if ~exist(OUT_ROOT,'dir'); mkdir(OUT_ROOT); end

for p=1:numel(ABP_FILES)
    PATH_ABP=ABP_FILES{p}; PATH_NUM=NUM_FILES{p};
    [~,abp_base,~]=fileparts(PATH_ABP);
    pat_label=sprintf('patient_%02d_%s',p,abp_base);
    OUT_DIR_P5=fullfile(OUT_ROOT,pat_label);
    if ~exist(OUT_DIR_P5,'dir'); mkdir(OUT_DIR_P5); end
    fprintf('\n================ %s =================\n',pat_label);

    raw=readmatrix(PATH_ABP,'FileType','text');
    t=raw(:,1); abp=raw(:,2); t=t(:); abp=abp(:);
    mask12=(t<=TMAX_SEC); time12=t(mask12); ABP12=abp(mask12);

    t_on=wabp(ABP12);
    feat=abpfeature(ABP12,t_on);
    beatq=jSQI(feat,t_on(1:end-1),ABP12);
    time=time12; ABP=ABP12;
    mat_out_q5=fullfile(OUT_DIR_P5,'first12h_Q5.mat');
    save(mat_out_q5,'time','ABP','t_on','feat','beatq');

    [co_uncal,to_min,~,fea]=estimateCO_v2(mat_out_q5,14,1);
    to_hr=to_min/60;
    Psys=fea(:,2); Pdias=fea(:,4); PP=fea(:,5); MAP=fea(:,6);
    Period=fea(:,7); HR=60*Fs./Period;

    tbl=readtable(PATH_NUM,'FileType','text','HeaderLines',2);
    time_num_all=tbl{:,1}; CO_TD_all=tbl{:,end};
    nz_idx_all=find(CO_TD_all~=0 & ~isnan(CO_TD_all));
    td_time_sec_all=time_num_all(nz_idx_all);
    td_CO_Lmin_all=CO_TD_all(nz_idx_all);
    mask12_TD=td_time_sec_all<=TMAX_SEC;
    td_time_sec_12h=td_time_sec_all(mask12_TD);
    td_CO_Lmin_12h=td_CO_Lmin_all(mask12_TD);
    td_time_hr_12h=td_time_sec_12h/3600;

    % --- C2 ---
    TD0=td_CO_Lmin_all(1); t0_sec=td_time_sec_all(1);
    [~,i_near]=min(abs(to_min-(t0_sec/60)));
    est0=co_uncal(i_near); k_C2=TD0/est0;

    % --- C1 (safe global LSQ) ---
    nTD=numel(td_time_sec_all); x_td=nan(nTD,1);
    for ii=1:nTD
        td_min=td_time_sec_all(ii)/60;
        [~,jj]=min(abs(to_min-td_min));
        if isempty(jj)||jj<1||jj>numel(co_uncal)
            x_td(ii)=NaN;
        else
            x_td(ii)=co_uncal(jj);
        end
    end
    valid_mask=~isnan(x_td)&isfinite(td_CO_Lmin_all);
    x_td_valid=x_td(valid_mask);
    td_CO_valid=td_CO_Lmin_all(valid_mask);
    if isempty(x_td_valid)
        warning('C1 calibration: no valid TD matches. Using k=1.');
        k_C1=1;
    else
        num_C1=sum(td_CO_valid.*x_td_valid);
        den_C1=sum(x_td_valid.^2);
        k_C1=num_C1/den_C1;
    end
    fprintf('Calibration factors: C2=%.3f, C1=%.3f (used %d valid TD points)\n', ...
        k_C2,k_C1,numel(x_td_valid));

    co_cal_C2=k_C2*co_uncal;
    co_cal_C1=k_C1*co_uncal;

    % --- Align time and CO arrays to same length before plotting ---
    L = min([length(to_hr), length(co_cal_C2), length(co_cal_C1)]);
    to_hr       = to_hr(1:L);
    co_cal_C2   = co_cal_C2(1:L);
    co_cal_C1   = co_cal_C1(1:L);

    % --- Plot CO (C1 vs C2) ---
    figure('Color','w','Name',sprintf('Parlikar - %s',pat_label));
    hold on;
    plot(to_hr,co_cal_C2,'Color',[0.5 0.5 0.5],'LineWidth',1.5,'DisplayName','C2 (1st TD)');
    plot(to_hr,co_cal_C1,'b--','LineWidth',1.2,'DisplayName','C1 (Global)');
    if ~isempty(td_time_hr_12h)
        stem(td_time_hr_12h,td_CO_Lmin_12h,'k','filled','DisplayName','Thermodilution');
    end
    ylabel('CO [L/min]'); xlabel('Time [hours]');
    title(sprintf('Parlikar (#14) – %s',pat_label),'Interpreter','none');
    legend('Location','best'); grid on; xlim([0 TMAX_HR]);
    saveas(gcf,fullfile(OUT_DIR_P5,'Problem5_Parlikar_C1C2.png')); close(gcf);
end
fprintf('\n=== Problem 5 Completed with C1 and C2 Calibrations ===\n');


%% -------------------------- PROBLEM 6 ----------------------------------
fprintf('\n=== Starting Problem 6: TPR Parlikar-style plots ===\n');
OUT_Q6_ROOT=fullfile(OUT_ROOT,'Q6');
if ~exist(OUT_Q6_ROOT,'dir'); mkdir(OUT_Q6_ROOT); end
PATS_FOR_Q6=[1 3 4];

for idx=1:numel(PATS_FOR_Q6)
    p=PATS_FOR_Q6(idx);
    PATH_ABP=ABP_FILES{p}; PATH_NUM=NUM_FILES{p};
    [~,abp_base,~]=fileparts(PATH_ABP);
    pat_label=sprintf('patient_%02d_%s',p,abp_base);
    OUT_DIR_Q5=fullfile(OUT_ROOT,pat_label);
    OUT_DIR_Q6=fullfile(OUT_Q6_ROOT,pat_label);
    if ~exist(OUT_DIR_Q6,'dir'); mkdir(OUT_DIR_Q6); end
    fprintf('\n--- Generating TPR plot for %s ---\n',pat_label);

    mat_out_q5=fullfile(OUT_DIR_Q5,'first12h_Q5.mat');
    if ~isfile(mat_out_q5); error('Missing %s. Run Problem 5 first.',mat_out_q5); end
    [co_uncal,to_min,~,fea]=estimateCO_v2(mat_out_q5,14,1);
    to_sec=to_min*60; MAP=fea(:,6); Period=fea(:,7);
    HR=60*Fs./Period; %#ok<NASGU>

    tbl=readtable(PATH_NUM,'FileType','text','HeaderLines',2);
    time_num_all=tbl{:,1}; CO_TD_all=tbl{:,end};
    nz_idx_all=find(CO_TD_all~=0 & ~isnan(CO_TD_all));
    td_time_sec_all=time_num_all(nz_idx_all);
    td_CO_Lmin_all=CO_TD_all(nz_idx_all);
    if isempty(td_time_sec_all)
        k_cal=1;
    else
        TD0=td_CO_Lmin_all(1); t0_sec=td_time_sec_all(1);
        [~,i_near]=min(abs(to_min-(t0_sec/60)));
        if isempty(i_near)||i_near<1||i_near>numel(co_uncal)
            k_cal=1;
        else
            k_cal=TD0/co_uncal(i_near);
        end
    end
    co_cal=k_cal*co_uncal;
    L=min([length(to_sec),length(co_cal),length(MAP)]);
    to_sec=to_sec(1:L); co_cal=co_cal(1:L); MAP=MAP(1:L);
    CO_mlps=co_cal*1000/60; TPR=MAP./CO_mlps; TPR(~isfinite(TPR))=NaN;
    f=figure('Color','w','Name',sprintf('TPR_%s',pat_label));
    plot(to_sec,TPR,'b','LineWidth',1.2);
    ylabel('Resistance (mmHg/(mL/s))'); xlabel('Time (secs)');
    title(sprintf('Estimated TPR – %s',pat_label),'Interpreter','none');
    grid on; xlim([0 TMAX_SEC]);
    saveas(f,fullfile(OUT_DIR_Q6,sprintf('%s_TPR_only.png',pat_label))); close(f);
end
fprintf('\n=== Problem 6 (TPR-only) Completed ===\n');
