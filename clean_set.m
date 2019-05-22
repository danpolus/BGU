%
%Step 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

bad_chan_in_epoch_percent = 0.1; %percent of bad channels in epoch to consider it bad
epoch_noise_zScore_th = 7; %-> play with it .  channel zscore in epoch, to consider channel bad. Should be at least twice than avalanche detection zThresh
minimal_nof_bad_pnts_epoch_chan = 1; %minimal number of threshold crossed points in epoch in channel in order to reject it
%
bad_channel_zScore_th = 3; 
bad_channel_time_percent = 0.15; %part of bad data in channel to mark it bad 
minimal_interchannel_correlation = 0.6;
%
LOW_PASS_HZ = 45;
fp = 'D:\My Files\Work\BGU\datasets\Panas\';
CHANNEL_LOCATION_FILE_INTERPOLATE = 'D:\My Files\Work\BGU\scripts\Mental-Imagery\electrodes\chanlocs60.sfp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn, fp] = uigetfile([fp '*.set'], 'Select eeglab file');
EEG = pop_loadset([fp fn]);

%low-pass filter
EEG = pop_eegfiltnew(EEG, [], LOW_PASS_HZ);

%CLEAN NOISE
EEG.bad_channels = [];                            
EEG.reject.rejglobal = zeros(1,EEG.trials);
EEG.reject.rejglobalE = zeros(EEG.nbchan,EEG.trials);
EEG = pop_reref(EEG,[]);
EEG_epochsConcat = eeg2epoch(EEG);

%find bad chanels: badly correlated with other chanels
[EEG_clean,~,~] = clean_artifacts(EEG_epochsConcat, ...
                                'ChannelCriterion', minimal_interchannel_correlation,... 
                                'ChannelCriterionMaxBadTime', bad_channel_time_percent*2 ,...
                                'NoLocsChannelCriterion', 'off',...
                                'NoLocsChannelCriterionExcluded', 'off',...
                                'LineNoiseCriterion',  'off',...
                                'FlatlineCriterion', 'off',...
                                'BurstCriterion',    'off',... % epoch_noise_zScore_th,...
                                'BurstCriterionRefMaxBadChns', 'off',...
                                'BurstCriterionRefTolerances', 'off',...
                                'WindowCriterion',   'off',... %bad_chan_in_epoch_percent,...
                                'WindowCriterionTolerances', 'off',... %[-epoch_noise_zScore_th epoch_noise_zScore_th],...
                                'Highpass',          'off');
if isfield(EEG_clean.etc, 'clean_channel_mask')
    EEG.bad_channels = unique([EEG.bad_channels find(EEG_clean.etc.clean_channel_mask == 0)']);
end
% if isfield(EEG_clean.etc, 'clean_sample_mask')
%     for ep=1:EEG.trials
%         if sum(EEG_clean.etc.clean_sample_mask((ep-1)*EEG.pnts+1 : ep*EEG.pnts) == 0) > 0
%             EEG.reject.rejglobal(ep) = 1;
%             EEG.reject.rejglobalE(:,ep) = 1;
%         end
%     end
% end

%find bad chanels: high std
EEG.bad_channels = unique([EEG.bad_channels find_bad_cahnnels(EEG_epochsConcat.data, EEG_epochsConcat.srate, 'window_step_length_sec', EEG.pnts/EEG.srate,...
    'bad_time_percent', bad_channel_time_percent, 'max_amp_score_thresh', bad_channel_zScore_th, 'min_amp_score_thresh', 0)]);

%find bad chanels in epochs
% EEG_ = pop_eegfiltnew(EEG, 15,[]); %use high-passed data
bad_epoch_chan = find_bad_epoch_channels(EEG.data, 'minimal_nof_bad_pnts', minimal_nof_bad_pnts_epoch_chan, 'z_score_thresh', epoch_noise_zScore_th);
bad_epoch_chan(EEG.bad_channels,:) = 0;
EEG.reject.rejglobal(sum(bad_epoch_chan,1)/(EEG.nbchan-length(EEG.bad_channels)) > bad_chan_in_epoch_percent) = 1;
EEG.reject.rejglobalE(bad_epoch_chan > 0) = 1;
EEG.reject.rejglobalE(EEG.bad_channels,:) = 1;

%plot
channels_color = repmat({'b'},1,size(EEG.data,1));   
channels_color(EEG.bad_channels) = {'r'};
EEG.reject.rejmanual = EEG.reject.rejglobal;
EEG.reject.rejmanualE = EEG.reject.rejglobalE;
pop_eegplot(EEG, 1, 0, 0, [], 'srate',EEG.srate, 'winlength',10, 'color', channels_color, 'spacing',60, 'eloc_file',[]);
%

%remove bad epochs
EEG.reject_hstr.rejglobal = EEG.reject.rejglobal;
EEG.reject_hstr.rejglobalE = EEG.reject.rejglobalE;
EEG = pop_rejepoch(EEG, EEG.reject.rejglobal ,0);

% %for calculation others than Avalanche Detection, interpolate bad channels
% EEG = pop_select( EEG, 'nochannel', EEG.bad_channels);
% EEG = eeg_interp(EEG, readlocs(CHANNEL_LOCATION_FILE_INTERPOLATE));

EEG = pop_reref(EEG,[]);

EEG = pop_saveset(eeg_checkset(EEG), 'filename',[EEG.filepath '\' EEG.filename(1:end-4) '_FiltClean.set'], 'savemode','onefile');
