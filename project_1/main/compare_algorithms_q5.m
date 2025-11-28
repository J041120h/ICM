function compare_algorithms_q5(ABP_FILES, NUM_FILES, OUT_ROOT, alg_ids, use_existing_mat)
%COMPARE_ALGORITHMS_Q5  Quantitatively compare CO-from-ABP algorithms at TD times (stable + debug-safe).
%
% Usage:
%   compare_algorithms_q5(ABP_FILES, NUM_FILES, OUT_ROOT);
%   compare_algorithms_q5(ABP_FILES, NUM_FILES, OUT_ROOT, [14 5 2 7], true);
%
% Notes:
%   - Uses C2 single-point calibration (first nonzero TD) per patient/algorithm.
%   - Robust interpolation (handles NaNs, duplicates, non-monotone x), with fallbacks.
%   - **Uniform row template** prevents "Subscripted assignment between dissimilar structures."
%
% Outputs:
%   <OUT_ROOT>/summary_q5/
%     - metrics_<patient>.csv
%     - metrics_all_patients.csv
%     - metrics_summary_macro_mean_sd.csv

if nargin < 4 || isempty(alg_ids),          alg_ids = [14 5 2 7]; end
if nargin < 5 || isempty(use_existing_mat), use_existing_mat = true; end

assert(iscellstr(ABP_FILES) && iscellstr(NUM_FILES), 'ABP_FILES and NUM_FILES must be cellstr.');
assert(numel(ABP_FILES) == numel(NUM_FILES), 'ABP_FILES and NUM_FILES sizes differ.');

if ~exist(OUT_ROOT,'dir'), mkdir(OUT_ROOT); end
SUM_DIR = fullfile(OUT_ROOT, 'summary_q5');
if ~exist(SUM_DIR,'dir'), mkdir(SUM_DIR); end

Fs       = 125;
TMAX_HR  = 12;
TMAX_SEC = TMAX_HR * 3600;

alg_name = containers.Map( ...
  {2, 5, 7, 14}, ...
  {'Pulse Pressure (#2)', 'Liljestrand (#5)', 'Corrected Impedance (#7)', 'Parlikar (#14)'} ...
);

fprintf('\n[Q5] === Starting comparison across %d patients ===\n', numel(ABP_FILES));

% --------- define a UNIFORM row template (all fields present) ----------
rowTemplate = struct( ...
    'Patient',            "", ...
    'AlgID',              NaN, ...
    'Algorithm',          "", ...
    'N_TD',               NaN, ...
    'RMSE',               NaN, ...
    'MAE',                NaN, ...
    'Bias',               NaN, ...
    'SD',                 NaN, ...
    'LOA_Lower',          NaN, ...
    'LOA_Upper',          NaN, ...
    'Pearson_r',          NaN, ...
    'DirAgreePct',        NaN, ...
    'DirPairs',           NaN, ...
    'PctAbsErr_le_0_5',   NaN, ...
    'PctAbsErr_le_1_0',   NaN, ...
    'CCC',                NaN, ...
    'k_cal',              NaN, ...
    'InterpInfo',         "" ...
);

all_rows = rowTemplate([]);  % empty struct array with fixed schema

for p = 1:numel(ABP_FILES)
    PATH_ABP = ABP_FILES{p};
    PATH_NUM = NUM_FILES{p};

    [~, abp_base, ~] = fileparts(PATH_ABP);
    pat_label = sprintf('patient_%02d_%s', p, abp_base);
    PAT_DIR   = fullfile(OUT_ROOT, pat_label);
    if ~exist(PAT_DIR,'dir'), mkdir(PAT_DIR); end

    fprintf('\n[Q5][%s] Prepare 12h context & TD\n', pat_label);

    % ---- Ensure first-12h .mat ----
    mat_path = fullfile(PAT_DIR, 'first12h_Q5.mat');
    if ~(use_existing_mat && exist(mat_path,'file'))
        raw = readmatrix(PATH_ABP,'FileType','text');
        if size(raw,2) < 2, error('ABP file must have at least 2 columns [time, ABP].'); end
        t   = raw(:,1);   abp = raw(:,2);
        mask12 = t <= TMAX_SEC;
        time = t(mask12);
        ABP  = abp(mask12);
        t_on  = wabp(ABP);
        feat  = abpfeature(ABP, t_on);
        beatq = jSQI(feat, t_on(1:end-1), ABP);
        if size(feat,1) - 1 ~= size(beatq,1)
            M = min(size(beatq,1), size(feat,1)-1);
            beatq = beatq(1:M,:);
        end
        save(mat_path,'time','ABP','t_on','feat','beatq');
    end

    % ---- Thermodilution (first 12h) ----
    tbl = readtable(PATH_NUM,'FileType','text','HeaderLines',2);
    time_num_all = tbl{:,1};     % seconds
    CO_TD_all    = tbl{:,end};   % L/min
    nz_idx_all      = find(CO_TD_all ~= 0 & ~isnan(CO_TD_all));
    td_time_sec_all = time_num_all(nz_idx_all);
    td_CO_Lmin_all  = CO_TD_all(nz_idx_all);

    mask12_TD       = td_time_sec_all <= TMAX_SEC;
    td_time_sec_12h = td_time_sec_all(mask12_TD);
    td_CO_Lmin_12h  = td_CO_Lmin_all(mask12_TD);
    td_time_min_12h = td_time_sec_12h / 60;

    fprintf('[Q5][%s] TD points in first 12h: %d\n', pat_label, numel(td_time_min_12h));
    if isempty(td_time_min_12h)
        warning('[Q5][%s] No nonzero TD in first 12h. Skipping patient.', pat_label);
        continue;
    end

    rows = rowTemplate([]);  % per-patient accumulator

    % ---------- Run each algorithm ----------
    for a = 1:numel(alg_ids)
        alg = alg_ids(a);
        name = getOr(alg_name, alg, sprintf('Alg #%d',alg));
        fprintf('[Q5][%s] estimateCO_v2 alg=%d (%s)\n', pat_label, alg, name);

        % Run estimator
        try
            [co_uncal,to_min,~,~] = estimateCO_v2(mat_path, alg, 1);
        catch ME
            warn_and_dump(ME, PAT_DIR, pat_label, 'estimateCO_v2_fail', struct('alg',alg,'mat_path',mat_path));
            continue;
        end

        % C2 single-point calibration at first TD
        TD0    = td_CO_Lmin_all(1);
        t0_sec = td_time_sec_all(1);
        [~,i_near] = min(abs(to_min - (t0_sec/60)));
        if isempty(i_near) || i_near<1 || i_near>numel(co_uncal) || ~isfinite(co_uncal(i_near))
            k_cal = 1;
        else
            k_cal = TD0 / co_uncal(i_near);
        end
        co_cal = k_cal * co_uncal;

        % Interpolate calibrated estimate at TD times (robust)
        x = to_min(:);
        y = co_cal(:);
        tq = td_time_min_12h(:);
        [y_at_td, dbg] = safe_interp_at(x, y, tq, PAT_DIR, pat_label, alg, name);
        if isempty(y_at_td)
            fprintf('[Q5][%s][%s] Interp failed; skipping.\n', pat_label, name);
            continue;
        end

        % Metrics
        td = td_CO_Lmin_12h(:);
        if numel(td) ~= numel(y_at_td)
            L = min(numel(td), numel(y_at_td));
            td = td(1:L); y_at_td = y_at_td(1:L);
        end

        e  = y_at_td - td;

        row = rowTemplate;  % start from uniform template
        row.Patient            = string(pat_label);
        row.AlgID              = alg;
        row.Algorithm          = string(name);
        row.N_TD               = numel(td);
        row.RMSE               = sqrt(mean(e.^2));
        row.MAE                = mean(abs(e));
        row.Bias               = mean(e);
        row.SD                 = std(e);
        row.LOA_Lower          = row.Bias - 1.96*row.SD;
        row.LOA_Upper          = row.Bias + 1.96*row.SD;
        row.Pearson_r          = corr_safe(y_at_td, td);
        [row.DirAgreePct, row.DirPairs] = directional_agreement(y_at_td, td);
        row.PctAbsErr_le_0_5   = 100*mean(abs(e) <= 0.5);
        row.PctAbsErr_le_1_0   = 100*mean(abs(e) <= 1.0);
        row.CCC                = concordance_cc(y_at_td, td);
        row.k_cal              = k_cal;
        row.InterpInfo         = string(dbg);

        rows(end+1) = row; %#ok<AGROW>
    end

    % Per-patient CSV
    if ~isempty(rows)
        T = struct2table(rows);
        per_csv = fullfile(SUM_DIR, sprintf('metrics_%s.csv', pat_label));
        writetable(T, per_csv);
        fprintf('[Q5][%s] Saved metrics: %s\n', pat_label, per_csv);

        all_rows = [all_rows; rows(:)]; %#ok<AGROW>
    else
        fprintf('[Q5][%s] No rows produced for this patient.\n', pat_label);
    end
end

% Aggregate CSV + macro means
if ~isempty(all_rows)
    Tall = struct2table(all_rows);
    agg_csv = fullfile(SUM_DIR, 'metrics_all_patients.csv');
    writetable(Tall, agg_csv);
    fprintf('[Q5] Saved aggregate: %s\n', agg_csv);

    algs = unique(Tall.AlgID);
    Srows = rowTemplate([]); % temp holder with same schema, but we'll only populate summary fields into a struct array
    summary = struct( ...
        'Algorithm', "", 'AlgID', NaN, 'N_Patients', NaN, ...
        'RMSE_mean', NaN, 'RMSE_sd', NaN, ...
        'MAE_mean', NaN, 'MAE_sd', NaN, ...
        'Bias_mean', NaN, 'Bias_sd', NaN, ...
        'SD_mean', NaN, 'SD_sd', NaN, ...
        'r_mean', NaN, 'r_sd', NaN, ...
        'CCC_mean', NaN, 'CCC_sd', NaN, ...
        'LOAlo_mean', NaN, 'LOAlo_sd', NaN, ...
        'LOAhi_mean', NaN, 'LOAhi_sd', NaN, ...
        'DirAgree_mean', NaN, 'DirAgree_sd', NaN, ...
        'pct05_mean', NaN, 'pct05_sd', NaN, ...
        'pct10_mean', NaN, 'pct10_sd', NaN ...
    );
    Sacc = summary([]); %#ok<NASGU>
    Sacc = summary([]); % ensure empty with schema

    SS = summary([]); % collector for summaries
    for i = 1:numel(algs)
        A = Tall(Tall.AlgID==algs(i), :);
        S = summary;
        S.Algorithm   = string(A.Algorithm(1));
        S.AlgID       = algs(i);
        S.N_Patients  = numel(unique(A.Patient));
        [S.RMSE_mean, S.RMSE_sd] = mstats(A.RMSE);
        [S.MAE_mean,  S.MAE_sd ] = mstats(A.MAE);
        [S.Bias_mean, S.Bias_sd] = mstats(A.Bias);
        [S.SD_mean,   S.SD_sd  ] = mstats(A.SD);
        [S.r_mean,    S.r_sd   ] = mstats(A.Pearson_r);
        [S.CCC_mean,  S.CCC_sd ] = mstats(A.CCC);
        [S.LOAlo_mean,S.LOAlo_sd] = mstats(A.LOA_Lower);
        [S.LOAhi_mean,S.LOAhi_sd] = mstats(A.LOA_Upper);
        [S.DirAgree_mean, S.DirAgree_sd] = mstats(A.DirAgreePct);
        [S.pct05_mean, S.pct05_sd] = mstats(A.PctAbsErr_le_0_5);
        [S.pct10_mean, S.pct10_sd] = mstats(A.PctAbsErr_le_1_0);
        SS(end+1) = S; %#ok<AGROW>
    end
    TS = struct2table(SS);
    sum_csv = fullfile(SUM_DIR, 'metrics_summary_macro_mean_sd.csv');
    writetable(TS, sum_csv);
    fprintf('[Q5] Saved macro summary: %s\n', sum_csv);
else
    warning('[Q5] No aggregate rows produced.');
end

fprintf('[Q5] === Done. ===\n');

end % ===== main =====


% ========================== Helper Functions ===========================

function [yq, dbg] = safe_interp_at(x, y, tq, OUT_DIR, pat_label, alg, alg_name)
dbg = "ok"; yq = [];
x = x(:); y = y(:); tq = tq(:);

% Equalize lengths
Lxy = min(numel(x), numel(y));
x = x(1:Lxy); y = y(1:Lxy);

% Finite only
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

if numel(x) < 2, return; end

% Deduplicate & sort strictly increasing
[xu, ia, ic] = unique(x, 'stable');
if numel(xu) < numel(x)
    y = accumarray(ic, y, [], @mean);
    x = xu;
    dbg = dbg + "|dedup";
end
% Ensure strictly increasing; if not, sort
if any(diff(x) <= 0)
    [x, idx] = sort(x, 'ascend');
    y = y(idx);
    dbg = dbg + "|sorted";
end

% Clip queries to domain
tq_clip = tq;
tq_clip(tq_clip < x(1))   = x(1);
tq_clip(tq_clip > x(end)) = x(end);

% Try linear, then nearest, then explicit nearest
try
    yq = interp1(x, y, tq_clip, 'linear');
    bad = ~isfinite(yq);
    if any(bad)
        yq(bad) = interp1(x, y, tq_clip(bad), 'nearest', 'extrap');
        dbg = dbg + "|fillNearest";
    end
catch
    try
        yq = interp1(x, y, tq_clip, 'nearest', 'extrap');
        dbg = dbg + "|fallbackNearest";
    catch
        % explicit nearest
        yq = zeros(size(tq_clip));
        for i = 1:numel(tq_clip)
            [~,k] = min(abs(x - tq_clip(i)));
            yq(i) = y(k);
        end
        dbg = dbg + "|explicitNearest";
    end
end
end

function warn_and_dump(ME, OUT_DIR, pat_label, tag, extra)
fprintf('[Q5][%s] %s: %s\n', pat_label, tag, ME.message);
dump_path = fullfile(OUT_DIR, sprintf('DEBUG_%s_%s.mat', pat_label, tag));
try, save(dump_path, 'ME', 'extra'); catch, end
end

function r = corr_safe(x, y)
x = x(:); y = y(:);
if numel(x) < 2 || numel(y) < 2 || all(~isfinite(x)) || all(~isfinite(y))
    r = NaN; return;
end
try
    r = corr(x, y, 'rows', 'pairwise');
catch
    r = NaN;
end
end

function [pct, n_pairs] = directional_agreement(yhat, y)
d1 = diff(y); d2 = diff(yhat);
keep = isfinite(d1) & isfinite(d2) & (d1 ~= 0);
n_pairs = sum(keep);
if n_pairs == 0
    pct = NaN;
else
    pct = 100 * mean(sign(d1(keep)) == sign(d2(keep)));
end
end

function ccc = concordance_cc(x, y)
x = x(:); y = y(:);
mu_x = mean(x); mu_y = mean(y);
s_x2 = var(x,1); s_y2 = var(y,1);
s_xy = mean((x-mu_x).*(y-mu_y));
ccc  = (2*s_xy) / (s_x2 + s_y2 + (mu_x - mu_y)^2);
end

function [m,s] = mstats(v)
v = v(:);
m = mean(v,'omitnan');
s = std(v,'omitnan');
end

function val = getOr(m, key, defaultVal)
if isKey(m, key), val = m(key); else, val = defaultVal; end
end
