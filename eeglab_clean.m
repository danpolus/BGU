%
%Step 1 : eeglab_clean - clean subject data with eeglab
%
%  inputs:
% fp
% fn - mat file name
% interpolateBadChannelsFlg
% saveFlg
% plotFlg
%
%  outputs:
% EEGSets - cleaned data separated into words and conditions, as EEGlab sets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EEGSets = eeglab_clean(fp, fn, interpolateBadChannelsFlg, saveFlg, plotFlg)

CHANNEL_LOCATION_FILE = 'D:\My Files\Work\BGU\scripts\Mental-Imagery\electrodes\chanlocs64.sfp';
CHANNEL_LOCATION_FILE_INTERPOLATE = 'D:\My Files\Work\BGU\scripts\Mental-Imagery\electrodes\chanlocs60.sfp';
%
dataset_name = 'eeg_data_wrt_task_rep_no_eog_256Hz_';
conditions = {'end_trial', 'last_beep'};
%
nof_chan = 64;
EOG_CHANNELS = [1 10 33 64];
fs = 256;
%
LOW_PASS_HZ = 45;
%
bad_chan_in_epoch_percent = 0.1; %percent of bad channels in epoch to consider it bad
epoch_noise_zScore_th = 7; %-> play with it .  channel zscore in epoch, to consider channel bad. Should be at least twice than avalanche detection zThresh
minimal_nof_bad_pnts_epoch_chan = 1; %minimal number of threshold crossed points in epoch in channel in order to reject it
%
bad_channel_zScore_th = 3;
bad_channel_time_percent = 0.15; %part of bad data in channel to mark it bad
minimal_interchannel_correlation = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%

EEGSets = [];

k = strfind(fp,'\');
FileInfo.base_fp = fp(1:k(end-1));
FileInfo.orig_fn = fn(1:end-4);
if contains(fp,'Long_words')
    FileInfo.scenario = '1LongWords';
elseif contains(fp,'Short_long_words')
    FileInfo.scenario = '2ShortLongWords';
elseif contains(fp,'Short_words')
    FileInfo.scenario = '3ShortWords';
elseif contains(fp,'Vowels')
    FileInfo.scenario = '4Vowels';
end
FileInfo.subj_id = fn(5:6);

if saveFlg
    output_fp = [FileInfo.base_fp '1 eeglab\'];
    mkdir(output_fp);
end

eeg_data_wrt = load(fullfile(fp,fn));

for iCond = 1:length(conditions)
    
    if iCond == 1
        FileInfo.condition = '1rest';
    elseif iCond == 2
        FileInfo.condition = '2imagine';
    else
        error('conditions');
    end    
    dataset = eval(['eeg_data_wrt.' dataset_name conditions{iCond}]);
    
    for iWord = 1:size(dataset,1)
        FileInfo.word_num = num2str(iWord);
        dataset_mat = cell2mat(dataset(iWord,:));
        
        %check correct nof channels
        if size(dataset_mat,1) > nof_chan
            warning(['Subject has more than' num2str(nof_chan) ' channels']);
            dataset_mat = dataset_mat(1:nof_chan,:);
        elseif size(dataset_mat,1) < nof_chan
            error(['Subject has LESS than' num2str(nof_chan) ' channels!']);
        end
        
        %create eeglab dataset
        EEG = pop_importdata('setname',[FileInfo.scenario  ' ' conditions{iCond} ' word' num2str(iWord) ' - ' FileInfo.orig_fn], 'data',dataset_mat, 'dataformat','array', 'chanlocs',CHANNEL_LOCATION_FILE, 'nbchan',size(dataset_mat,1), 'pnts',size(dataset{1,1},2), 'srate',fs);
        
        EEG.FileInfo = FileInfo;
        
        %remove eog channels
        EEG = pop_select( EEG,'nochannel', EOG_CHANNELS);
        
%         %save berfore cleaning
%         if saveFlg
%             output_fn = [output_fp FileInfo.orig_fn '_' conditions{iCond} '_word' num2str(iWord) '_dirty'];
%             EEG = pop_saveset(eeg_checkset(EEG), 'filename',[output_fn '.set'], 'savemode','onefile');
%             csvwrite([output_fn '_ch' num2str(EEG.nbchan) '_pt' num2str(EEG.pnts) '_ep' num2str(EEG.trials) '_fs' num2str(EEG.srate) '.csv'], reshape(EEG.data,size(EEG.data,1),[]));
%         end
        
        %low-pass filter
        EEG = pop_eegfiltnew(EEG, [], LOW_PASS_HZ);
        
        %CLEAN NOISE
        EEG.bad_channels = [];
        EEG.reject.rejglobal = zeros(1,EEG.trials);
        EEG.reject.rejglobalE = zeros(EEG.nbchan,EEG.trials);
        EEG = pop_reref(EEG,[]);
        EEG_epochsConcat = eeg2epoch(EEG);
        
        %find bad chanels: badly correlated with other channels
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
        if plotFlg
            channels_color = repmat({'b'},1,size(EEG.data,1));
            channels_color(EEG.bad_channels) = {'r'};
            EEG.reject.rejmanual = EEG.reject.rejglobal;
            EEG.reject.rejmanualE = EEG.reject.rejglobalE;
            pop_eegplot(EEG, 1, 0, 0, [], 'srate',EEG.srate, 'winlength',10, 'color', channels_color, 'spacing',60, 'eloc_file',[]);
        end
        
        %remove bad epochs
        EEG.reject_hstr.rejglobal = EEG.reject.rejglobal;
        EEG.reject_hstr.rejglobalE = EEG.reject.rejglobalE;
        EEG = pop_rejepoch(EEG, EEG.reject.rejglobal ,0);
        
        %for calculation others than Avalanche Detection, interpolate bad channels
        if interpolateBadChannelsFlg && ~isempty(EEG.bad_channels)
            EEG = pop_select( EEG, 'nochannel', EEG.bad_channels);
            EEG = eeg_interp(EEG, readlocs(CHANNEL_LOCATION_FILE_INTERPOLATE));
        end
        
        EEG = pop_reref(EEG,[]);
        
        EEGSets = [EEGSets eeg_checkset(EEG)];

        %save eeglab set and csv file
        if saveFlg
            output_fn = [output_fp FileInfo.orig_fn '_' conditions{iCond} '_word' num2str(iWord)];
            EEG = pop_saveset(eeg_checkset(EEG), 'filename',[output_fn '.set'], 'savemode','onefile');
            csvwrite([output_fn '_ch' num2str(EEG.nbchan) '_pt' num2str(EEG.pnts) '_ep' num2str(EEG.trials) '_fs' num2str(EEG.srate) '.csv'], reshape(EEG.data,size(EEG.data,1),[]));
        end

    end
    
end
