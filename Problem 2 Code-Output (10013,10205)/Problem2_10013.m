% Project 1, Problem 2 (a) Code
% Reproduce traces similar to Figure 1 in paper for PATIENT 10013
% ONLY produce plots for 11h
clear; clc;

% --------------------------DEFINING VARIABLES------------------------
% define paths for inputs/outputs. EDIT before running
PATH_ABP = ['C:\Users\danie\OneDrive\Desktop\ICM\s10013\s10013\' ...
    's10013-2564-11-01-23-37_ABP.txt'];   % change for path of ABP file
OUT_DIR  = 'C:\Users\danie\OneDrive\Desktop\ICM\Project1_Output'; % path for output

% create output folder if it does not exist
if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

% define time anchors from the ABP records
% as we plot the first 20 pulses starting at 10 h and 11 h, we find
% timepoint for 10hours and 11hours in the ABP file.
T10 = 10*3600;   % 10 hours
T11 = 11*3600;   % 11 hours

% -------------------------------READ DATA--------------------------------
% read the ABP file: col 1 is time in seconds, col 2 is ABP in mmHg
raw = readmatrix(PATH_ABP,'FileType','text');
t   = raw(:,1);
abp = raw(:,2);
t   = t(:); % enforce vectors to be column vectors
abp = abp(:);

% ---------------------------BEAT ONSET DETECTION-------------------------
% obtain onset time for each beat using helper function
OnsetTimes = wabp(abp);

% ------------------- SELECT 21 ONSETS (FOR 20 BEATS) ---------------------
% General Algo:
%   1) Find the first sample index in the time vector at or after the anchor.
%   2) Find the first detected onset at or after that sample index.
%   3) Take 21 consecutive onset indices to ensure exactly 20 beats.

% for 10 h
start10_idx = find(t >= T10, 1, 'first');
k10 = find(OnsetTimes >= start10_idx, 1, 'first');
onsets10 = OnsetTimes(k10:(k10+20));   % 21 onsets = 20 beats

% For 11 h
start11_idx = find(t >= T11, 1, 'first');
k11 = find(OnsetTimes >= start11_idx, 1, 'first');
onsets11 = OnsetTimes(k11:(k11+20));   % 21 onsets = 20 beats

% --------------- EXTRACTING FEATURES -----------------------------------
% from abpfeature.m,
%   col  9: End of systole time by 0.3*sqrt(RR) method
%   col 11: End of systole time by lowest negative slope method
feat10 = abpfeature(abp, onsets10); 
feat11 = abpfeature(abp, onsets11);

% ----------------------- CALCULATE SIGNAL QUALITY ------------------------
% calculate ABP waveform signal quality index using helper function
[BeatQ10, ~] = jSQI(feat10, onsets10(1:end-1), abp);
[BeatQ11, ~] = jSQI(feat11, onsets11(1:end-1), abp);

% ------------------------------ PLOT (11 h) -------------------------------
% the logic is the exact same as above.
i0 = onsets11(1);
i1 = onsets11(end);
tt = t(i0:i1);
yy = abp(i0:i1);

on20   = onsets11(1:end-1);
eos03  = feat11(:,9);
eosMin = feat11(:,11);

f2 = figure('Color','w');
plot(tt, yy, 'LineWidth', 1); hold on;
plot(t(on20), abp(on20), '*', 'MarkerSize', 5);
plot(t(eos03), abp(eos03), 'xr', 'MarkerSize', 7, 'LineWidth', 1.2);
plot(t(eosMin), abp(eosMin), 'o', 'MarkerSize', 5);
xlabel('Time (s)');
ylabel('ABP (mmHg)');
title('Patient #10013, Starting at 11 h');
box on; grid on;

print(fullfile(OUT_DIR,'Problem2_10013_11h.png'),'-dpng','-r300');
close(gcf);