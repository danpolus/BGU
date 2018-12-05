%
%Step 1
%
% Instructions For Manual Rejection:
% ----------------------------------
% 1. If nof channels > 10%, reject the whole epoch
% 2. Absolute threshlolds should be the same for single subject scenario
% 3. Use manual channel epoch rejection for final adjustment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

signal_uV_th = 50; %absolute threshlolds should be the same for single subject scenario
spacing_uV = 100;
%
bad_channel_zScore_th = 3;
epoch_noise_zScore_th = 6;
%
bad_time_percent = 0.3; %part of bad data in channel to mark it bad 
LOW_PASS_HZ = 45;
fp = 'D:\My Files\Work\BGU\datasets\Panas\';
CHANNEL_LOCATION_FILE_INTERPOLATE = 'D:\My Files\Work\BGU\scripts\Mental Imagery\electrodes\chanlocs60.sfp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn, fp] = uigetfile([fp '*.set'], 'Select eeglab file');
EEG = pop_loadset([fp fn]);

%low-pass filter
EEG = pop_eegfiltnew(EEG, [], LOW_PASS_HZ);

EEG = pop_saveset(eeg_checkset(EEG), 'filename',[EEG.filepath '\' EEG.filename(1:end-4) '_FiltClean.set'], 'savemode','onefile');

EEG.bad_channels = [];

%auto bad channels
EEG_epochsConcat = pop_reref(EEG,[]);
EEG_epochsConcat.data = reshape(EEG_epochsConcat.data, size(EEG_epochsConcat.data,1), []);
EEG_epochsConcat.trials = 1;  EEG_epochsConcat.epoch = []; EEG_epochsConcat.pnts = size(EEG_epochsConcat.data,2);
auto_bad_channels = find_bad_cahnnels(EEG_epochsConcat.data, EEG_epochsConcat.srate, 'window_step_length_sec', EEG.pnts/EEG.srate, 'bad_time_percent', bad_time_percent/2, 'max_amp_score_thresh', bad_channel_zScore_th, 'min_amp_score_thresh', 1/bad_channel_zScore_th);
disp([EEG.setname '    suspected bad channels: ' num2str(auto_bad_channels) '?']);

% %manual bad channels
% csc_eeg_plotter(EEG_epochsConcat);
% manual_bad_channels = cell2mat(hidden_channels(:,2))';
% EEG.bad_channels = [EEG.bad_channels manual_bad_channels];
% %EEG = pop_select( EEG,'nochannel', manual_bad_channels);
% EEG = pop_saveset(eeg_checkset(EEG), 'filename',[EEG.filepath '\' EEG.filename], 'savemode','onefile');


%BAD EPOCHS AND CHANNELS
%reject extream uV values
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -signal_uV_th, signal_uV_th, EEG.xmin, EEG.xmax, 0, 0);
% %Reject abnormal trend
% EEG = pop_rejtrend(EEG,1,1:EEG.nbchan ,150,50,0.3,2,0);
% %Reject by epoch and channel z-score
% EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, epoch_noise_zScore_th, bad_channel_zScore_th, 0, 0, 0, [], 0);
%Reject channels kurtosis
% EEG = pop_rejkurt(EEG, 1, 1:EEG.nbchan, epoch_noise_zScore_th, bad_channel_zScore_th, 0, 0, 0, [], 0);
EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
%Update automatic rejection manualy
EEG.reject.rejmanual = EEG.reject.rejglobal;
EEG.reject.rejmanualE = EEG.reject.rejglobalE;
ALLEEG = EEG; CURRENTSET = 1;
pop_eegplot(EEG, 1, 0, 0, [], 'srate',EEG.srate, 'winlength',10, 'spacing', spacing_uV, 'eloc_file', []);
return;
%%%%
%Manual channel epoch rejection. ENTER ALL HERE:

ALLEEG.reject.rejglobalE(:,[16 78 97]) = 0;

bad_channels_manual = [];
EEG.bad_channels = [EEG.bad_channels bad_channels_manual];
%%%%
EEG.reject_hstr.rejmanual = ALLEEG.reject.rejmanual;
EEG.reject_hstr.rejglobalE = ALLEEG.reject.rejglobalE;
EEG = pop_rejepoch(EEG, EEG.reject.rejmanual ,0);
%%%
bad_epoch_chan = EEG.reject_hstr.rejglobalE(:,~EEG.reject_hstr.rejmanual);
% bad_channels_auto = find(mean(bad_epoch_chan,2) > bad_time_percent);
% EEG.bad_channels = [EEG.bad_channels bad_channels_auto];
EEG = pop_saveset(eeg_checkset(EEG), 'filename',[EEG.filepath '\' EEG.filename], 'savemode','onefile');


%BAD CHANNELS SPECTRA
EEG_epochsConcat = pop_reref(EEG,[]);
for iEpoch=1:EEG_epochsConcat.trials
    EEG_epochsConcat.data(bad_epoch_chan(:,iEpoch),:,iEpoch) = 0;
end
EEG_epochsConcat.data = reshape(EEG_epochsConcat.data, size(EEG_epochsConcat.data,1), []);
EEG_epochsConcat.trials = 1;  EEG_epochsConcat.epoch = []; EEG_epochsConcat.pnts = size(EEG_epochsConcat.data,2);
bad_channels_spectra = channel_map_topoplot(EEG_epochsConcat, [], true);
return
EEG.bad_channels = [EEG.bad_channels bad_channels_spectra];

%INTERPOLATE BAD CHANNELS
EEG.bad_channels = unique(EEG.bad_channels);
% EEG = pop_select( EEG, 'nochannel', EEG.bad_channels);
% EEG = eeg_interp(EEG, readlocs(CHANNEL_LOCATION_FILE_INTERPOLATE));
EEG = pop_saveset(eeg_checkset(EEG), 'filename',[EEG.filepath '\' EEG.filename], 'savemode','onefile');
% EEG = pop_saveset(eeg_checkset(EEG), 'savemode','resave');
