%
% fast *.set files evaluation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

LOW_PASS_HZ = 45;
fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.set'], 'Select eeglab files', 'MultiSelect','on');

fig_spect = figure();
fig_spect_lp = figure('Name','low pass 50Hz');
for iFile = 1 : length(files)
    
    EEG = pop_loadset([fp files{iFile}]);
    
    %spectra
    figure(fig_spect);
    subplot(ceil(length(files)/2),2,iFile); spect = pop_spectopo(EEG, 1, [0 floor(EEG.xmax*1000)], 'EEG' , 'percent', 15, 'freqrange',[0 floor(EEG.srate/2)],'electrodes','off');title(EEG.setname);
    %     figure;plot(spect(32:33,:)');
    
    %low-pass filter
    EEG = pop_eegfiltnew(EEG, [], LOW_PASS_HZ);
    figure(fig_spect_lp);
    subplot(ceil(length(files)/2),2,iFile); spect = pop_spectopo(EEG, 1, [0 floor(EEG.xmax*1000)], 'EEG' , 'percent', 15, 'freqrange',[0 floor(EEG.srate/2)],'electrodes','off');title(EEG.setname);
    
    % Automatic bad channels
    EEG_epochsConcat = pop_reref(EEG,[]);
    EEG_epochsConcat.data = reshape(EEG_epochsConcat.data, size(EEG_epochsConcat.data,1), []);
    EEG_epochsConcat.trials = 1;  EEG_epochsConcat.epoch = []; EEG_epochsConcat.pnts = size(EEG_epochsConcat.data,2);
    bad_channels_auto = find_bad_cahnnels(EEG_epochsConcat.data, EEG_epochsConcat.srate, 'window_step_length_sec', EEG.pnts/EEG.srate, 'bad_time_percent', 0.1, 'max_amp_score_thresh', 3, 'min_amp_score_thresh', 1/3);
    disp([EEG.setname '    suspected bad channels: ' num2str(bad_channels_auto) '?']);
    
end
