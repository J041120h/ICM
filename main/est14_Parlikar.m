
function x = est14_Parlikar(Period, MAP, Pdias, PP_in, varargin)
% EST14_PARLIKAR  Uncalibrated CO estimator from Parlikar et al. (2007).
%   x = est14_Parlikar(Period, MAP, Pdias, PP_in)
%
% Inputs (beat-wise vectors, length N):
%   Period : beat periods in SAMPLES (as provided by abpfeature)
%   MAP    : mean arterial pressure for each beat [mmHg]
%   Pdias  : diastolic arterial pressure for each beat [mmHg]
%   PP_in  : pulse pressure Psys - Pdias for each beat [mmHg] (optional; 
%            if unreliable, the method internally recomputes PP via alpha*(MAP - Pdias))
%
% Output:
%   x : uncalibrated CO (proportional units), length N-1 (aligned to beats 1..N-1)
%
% Method:
%   Implements the beat-averaged Windkessel relations from Parlikar et al. (2007):
%      (ΔP_n / T_n) + P̅_n / τ_n = PP_n / T_n   (Eq. 8)
%   Estimate 1/τ_n via local least squares over a sliding odd-length window (w),
%   then compute (uncalibrated) CO per Eq. (10) with C_n absorbed into calibration:
%      CO_n ∝ (ΔP_n / T_n) + P̅_n / τ_n
%
% References:
%   Parlikar TA, Heldt T, Ranade GV, Verghese GC. 
%   Model-Based Estimation of Cardiac Output and Total Peripheral Resistance.
%   Computers in Cardiology, 2007.
%
% Defaults / Notes:
%   - Fs is assumed 125 Hz (as in PhysioToolkit). Tn (s) = Period / Fs.
%   - ΔP_n is approximated using diastolic pressures at consecutive onsets:
%       ΔP_n ≈ DAP_{n+1} - DAP_{n}
%   - PP_n is estimated as alpha * (MAP_n - DAP_n) with alpha=2 (triangular pulse).
%     If PP_in is provided and non-empty, we blend: PP = 0.5*PP_in + 0.5*alpha*(MAP-Pdias).
%   - Output has length N-1 due to ΔP_n definition.
%
%   This function returns RELATIVE/UNCALIBRATED CO. Use calibration C2/C3 (Sun 2009/2005)
%   against thermodilution to obtain absolute L/min values.
%
%   Written for course project usage. 
%
%   Harry/Jiang + ChatGPT helper, 2025-11-05.

Fs = 125;                           % sampling rate [Hz]
alpha = 2;                          % PP ≈ alpha*(MAP - Pdias) per Parlikar Eq. (9)
w = 9;                              % odd window size for local LS fit of 1/tau
w = max(5, w + mod(w+1,2));         % ensure odd and at least 5

% enforce column vectors
Period = Period(:);
MAP    = MAP(:);
Pdias  = Pdias(:);
if nargin < 4 || isempty(PP_in)
    PP_in = [];
else
    PP_in = PP_in(:);
end

N = length(MAP);
if length(Period)~=N || length(Pdias)~=N
    error('Input vectors must have the same length.');
end
if N < 20
    error('Not enough beats to apply Parlikar estimator.');
end

% Compute per-beat quantities
Tn = Period / Fs;                   % seconds
DAP = Pdias;                        % use diastolic at onset as pressure foot
PP_alpha = alpha * (MAP - DAP);
if ~isempty(PP_in)
    PP = 0.5*PP_in + 0.5*PP_alpha;  % blend to be robust to reflections
else
    PP = PP_alpha;
end

% ΔP_n uses DAP at consecutive beats; drop last beat
dP = DAP(2:end) - DAP(1:end-1);
T  = Tn(1:end-1);
P  = MAP(1:end-1);
PPc = PP(1:end-1);

% Target for local LS: y = PP/T - dP/T  and solve y ≈ (1/τ)*P over a window
y = (PPc ./ T) - (dP ./ T);

% Preallocate inv_tau (aligned to beats 1..N-1)
inv_tau = nan(size(P));

halfw = floor(w/2);
for n = 1:length(P)
    i0 = max(1, n - halfw);
    i1 = min(length(P), n + halfw);
    Pwin = P(i0:i1);
    ywin = y(i0:i1);
    % Least squares for scalar slope (1/τ): minimize || ywin - slope*Pwin ||
    denom = sum(Pwin.^2);
    if denom > 0
        inv_tau(n) = sum(Pwin .* ywin) / denom;
    else
        inv_tau(n) = NaN;
    end
end

% Compute uncalibrated CO per Eq. (10) with C absorbed by later calibration
x = (dP ./ T) + P .* inv_tau;

% Clean obvious NaNs/Infs (edges)
bad = ~isfinite(x);
if any(bad)
    x(bad) = interp1(find(~bad), x(~bad), find(bad), 'linear', 'extrap');
end
end
