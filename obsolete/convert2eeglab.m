%
%Step 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

CHANNEL_LOCATION_FILE = 'D:\My Files\Work\BGU\scripts\Mental-Imagery\electrodes\chanlocs64.sfp';
fp = 'D:\My Files\Work\BGU\datasets\Panas\';
dataset_name = 'eeg_data_wrt_task_rep_no_eog_256Hz_';
conditions = {'end_trial', 'last_beep'};

nof_chan = 64;
EOG_CHANNELS = [1 10 33 64];
fs = 256;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select data files', 'MultiSelect','on');

k = strfind(fp,'\');
scenario = fp(k(end-1)+1:k(end)-1);

for iFile = 1:length(files)
    eeg_data_wrt = load(fullfile(fp,files{iFile}));
    for iCond = 1:length(conditions)
        dataset = eval(['eeg_data_wrt.' dataset_name conditions{iCond}]);
        for iWord = 1:size(dataset,1)
            dataset_mat = cell2mat(dataset(iWord,:));
            
            %check correct nuf channels
            if size(dataset_mat,1) > nof_chan
                warning(['Subject has more than' num2str(nof_chan) ' channels']);
                dataset_mat = dataset_mat(1:nof_chan,:);
            elseif size(dataset_mat,1) < nof_chan
                error(['Subject has LESS than' num2str(nof_chan) ' channels!']);
            end
            
            %create and save eeglab dataset
            EEG = pop_importdata('setname',[scenario ' ' conditions{iCond} ' word' num2str(iWord) ' - ' files{iFile}(1:end-4)], 'data','dataset_mat', 'dataformat','array', 'chanlocs',CHANNEL_LOCATION_FILE, 'nbchan',size(dataset_mat,1), 'pnts',size(dataset{1,1},2), 'srate',fs);
            EEG = pop_select( EEG,'nochannel', EOG_CHANNELS); %remove eog channels
            EEG = pop_saveset(eeg_checkset(EEG), 'filename',[fp(1:k(end-1)) 'eeglab\' files{iFile}(1:end-4) '_' conditions{iCond} '_word' num2str(iWord) '.set'], 'savemode','onefile');
        end
    end
end

% %DEBUG - time and frequency domain
% figure;plot(dataset{1,6}(3,:));
% chan_spec = abs(fft(dataset_mat(3,:)));
% f = 0:fs/length(chan_spec):fs/2;
% figure;plot(f,chan_spec(1:length(f)));
