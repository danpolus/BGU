%
%Step 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';
[files, fp] = uigetfile([fp '*.set'], 'Select all subject scenario files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end
EEG = pop_loadset([fp files{1}]);

%%%%%%%%%%%%%%%%%
global param

%Event Size - range of avalanche sizes
param.ES.min = 1;
param.ES.max = (floor(EEG.nbchan*2/100))*100;
param.ES.edges = unique(ceil(logspace(0,log10(param.ES.max),25))); % log spaced bins
%params
param.MaxMin = 'maxmin'; % one of 3 options: 'maxmin', 'max', 'min'. default: 'maxmin'
param.Fs = EEG.srate;
% times for raster plot
param.t1 = 4; % sec
param.t2 = 8; % sec
%more input params
max_tau_sec = 0.04; %40msec
param.tau_vec = 1/EEG.srate:1/EEG.srate:max_tau_sec; % array with delta t values for avalanche analysis
param.optimal_alpha = -1.5;
param.optimal_sigma = 1;

std_TH = 3.0; %set ZsubjTh3 to appears in folder name
% Zfile means that mean and std calculated per file
% Zsubj means that mean and std calculated over all subject files for the same scenario

zero_bad_epoch_channels_flg = true; %true -> BEfile_BC####, false ->BEnone_BCnone set to appear in folder name
zero_bad_chanels_flg = true ;% false -> BCnone set to appear in folder name
zero_subj_scenario_bad_channels_flg = false; %true -> BCsubj, false -> BCfile set to appear in folder name

%%%%%%%%%%%%%%%%%%%%

%calculate mean and std of the whole scenario for each channel
EEG_data_reshaped = [];
subj_scenario_bad_channels = [];
for iFile = 1:length(files)
    EEG = pop_loadset([fp files{iFile}]);
    bad_epoch_chan = EEG.reject_hstr.rejglobalE(:,~EEG.reject_hstr.rejmanual);
    for iEpoch=1:EEG.trials
        EEG.data(bad_epoch_chan(:,iEpoch),:,iEpoch) = NaN; %set bad epoch channels to NaN
    end
    EEG_data_reshaped = [EEG_data_reshaped reshape(EEG.data,size(EEG.data,1),[])];
    subj_scenario_bad_channels = [subj_scenario_bad_channels EEG.bad_channels];
end
subj_scenario_mean = mean(EEG_data_reshaped,2,'omitnan');
subj_scenario_std = std(EEG_data_reshaped,[],2,'omitnan');
subj_scenario_bad_channels = unique(subj_scenario_bad_channels);
clear('EEG_data_reshaped');

for iFile = 1:length(files)
    
    file_info = [];
    if contains(fp,'Long_words')
        file_info.scenario = '1LongWords';
    elseif contains(fp,'Short_long_words')
        file_info.scenario = '2ShortLongWords';
    elseif contains(fp,'Short_words')
        file_info.scenario = '3ShortWords';
    elseif contains(fp,'Vowels')
        file_info.scenario = '4Vowels';
    end
    file_info.subj_id = files{iFile}(5:6);
    if contains(files{iFile},'end_trial_word')
        file_info.condition = '1rest';
        file_info.word_num = files{iFile}(strfind(files{iFile},'end_trial_word') + length('end_trial_word'));
    elseif contains(files{iFile},'last_beep_word')
        file_info.condition = '2imagine';
        file_info.word_num = files{iFile}(strfind(files{iFile},'last_beep_word') + length('last_beep_word'));
    end
    
    EEG = pop_loadset([fp files{iFile}]);
    
    %zscore
    subj_scenario_mean_mat = repmat(subj_scenario_mean,[1,size(EEG.data,2),size(EEG.data,3)]);
    subj_scenario_std_mat = repmat(subj_scenario_std,[1,size(EEG.data,2),size(EEG.data,3)]);
    EEG.data = (EEG.data - subj_scenario_mean_mat)./ subj_scenario_std_mat;
    EEG.data(isinf(EEG.data)) = 0;
    
    %set bad data to 0
    if zero_bad_epoch_channels_flg
        bad_epoch_chan = EEG.reject_hstr.rejglobalE(:,~EEG.reject_hstr.rejmanual);
        if zero_bad_chanels_flg
            bad_epoch_chan(EEG.bad_channels,:) = 1;
            if zero_subj_scenario_bad_channels_flg
                bad_epoch_chan(subj_scenario_bad_channels,:) = 1;
            end
        end
        for iEpoch=1:EEG.trials
            EEG.data(bad_epoch_chan(:,iEpoch),:,iEpoch) = 0; %set bad epoch channels to 0
        end
    end
    
    %initialize results vector
    for iTau = 1:length(param.tau_vec)
        all_epochs(iTau).av_raster_epochs = [];
        all_epochs(iTau).av_raster = [];
        all_epochs(iTau).av_size_vec = [];
        all_epochs(iTau).av_dur_vec = [];
        all_epochs(iTau).sigma_vec = [];
    end
    %find avalanches
    for iEpoch = 1:EEG.trials   
        AvalancheResults = FindAvalanches(EEG.data(:,:,iEpoch),EEG.times,param.tau_vec,param.MaxMin,std_TH,1,0);
        for iTau = 1:length(param.tau_vec)
            all_epochs(iTau).av_raster_epochs(:,:,iEpoch) = AvalancheResults(iTau).av_raster;
            all_epochs(iTau).av_raster = [all_epochs(iTau).av_raster AvalancheResults(iTau).av_raster];
            all_epochs(iTau).av_size_vec = [all_epochs(iTau).av_size_vec AvalancheResults(iTau).av_size_vec];
            all_epochs(iTau).av_dur_vec = [all_epochs(iTau).av_dur_vec AvalancheResults(iTau).av_dur_vec];
            all_epochs(iTau).sigma_vec = [all_epochs(iTau).sigma_vec AvalancheResults(iTau).sigma_vec];
        end
    end    
    
    for iTau = 1:length(param.tau_vec)
        all_epochs(iTau).alpha = estimateParamML(param.ES.min,param.ES.max,'zeta',-param.optimal_alpha,all_epochs(iTau).av_size_vec);
        all_epochs(iTau).sigma = mean(all_epochs(iTau).sigma_vec);  
        % rate of events in each channel
        R(:,iTau) = sum(all_epochs(iTau).av_raster,2)/(size(all_epochs(iTau).av_raster,2)/EEG.srate); 
    end
    
    [~,tau_optimal_inx] = min(sqrt(([all_epochs.sigma] - param.optimal_sigma).^2 + (-[all_epochs.alpha] - param.optimal_alpha).^2));
    
    %save results
    save([fp files{iFile}(1:end-4) '_avalanches.mat'],'all_epochs','param','file_info');
    
    %plot
    xx = 1:300; yy = xx.^(param.optimal_alpha); 
    xSize = 7; ySize = 9.5; xLeft = (8.5 - xSize)/2; yTop = (11 - ySize)/2;%Graphic parameters
    for iTau = 1:length(param.tau_vec)
        
        if iTau~=tau_optimal_inx
            continue;
        end
      
        figure('Name',[EEG.setname '  \tau=' num2str(iTau) 'dt']);
        set(gcf,'Color','w');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
        set(gcf,'Position',[570 40 xSize*75 ySize*75])
        
        subplot(2,3,1);histogram(all_epochs(iTau).av_size_vec,max(all_epochs(iTau).av_size_vec));title(['Avalanche Size  \tau = ' num2str(param.tau_vec(iTau))]);
        subplot(2,3,4);histogram(all_epochs(iTau).av_dur_vec,max(all_epochs(iTau).av_dur_vec));title(['Avalanche Duration  \tau = ' num2str(param.tau_vec(iTau))]);
              
        subplot(2,3,2); RasterPlot(all_epochs(iTau).av_raster,param.Fs,param.tau_vec(iTau));title(['\tau = ' num2str(param.tau_vec(iTau)) ' - Raster']); drawnow
        subplot(2,3,3); RasterPlot(all_epochs(iTau).av_raster,param.Fs,param.tau_vec(iTau),[param.t1 param.t2]);title('zoom in'); drawnow
        subplot(2,3,5); stem(R(:,iTau));xlabel('Channel #');ylabel('Events/sec');title('Event rate in each channel');set(gca,'XLim',[1 size(all_epochs(iTau).av_raster,1)]);
        
        % calculate avalanche distribution
        [size_dist, ES_edges1] = NormalizedDist(all_epochs(iTau).av_size_vec,param.ES.edges); 
        subplot(2,3,6); loglog(ES_edges1,size_dist,'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
        set(gca,'XLim',[0 310]);
        xlabel('S');ylabel('P(S)');
        title(['\tau = ' num2str(1000*param.tau_vec(iTau)) ' ms' ', \alpha = ' num2str(all_epochs(iTau).alpha) ', \sigma = ' num2str(all_epochs(iTau).sigma)]);
 
    end   
    
    figure('Name',EEG.setname);scatter([all_epochs.sigma],-[all_epochs.alpha],'x');xlabel('\sigma');ylabel('\alpha');title('\alpha = F(\sigma)  labels: \tau values [sec]');
    hold on
    scatter([all_epochs(tau_optimal_inx).sigma],-[all_epochs(tau_optimal_inx).alpha],'dk','filled');
    text([all_epochs.sigma]+0.01,-[all_epochs.alpha],strsplit(num2str(param.tau_vec)));
    plot([min(param.optimal_sigma,min([all_epochs.sigma]))-0.1, max(param.optimal_sigma,max([all_epochs.sigma]))+0.1],[param.optimal_alpha,param.optimal_alpha],'--m',...
        [param.optimal_sigma,param.optimal_sigma],[min(param.optimal_alpha,min(-[all_epochs.alpha]))-0.1, max(param.optimal_alpha,max(-[all_epochs.alpha]))+0.1],'--m');
    hold off
    
    
%     %find avalanche distribution for different nof electrodes
%
%     electrodes{1} = 1:60; description{1} = 'all electrodes';
%     electrodes{2} = 1:31; description{2} = 'every other electrode';
%     electrodes{3} = 1:2:31; description{3} = 'every forth electrode';
%     electrodes{4} = [6 37 27 38 22 55 10 51 21]; description{4} = 'small central spot electrodes';
%     electrodes{5} = [electrodes{4} 2 33 1 59 28 56 23 52 17 50 11 42 12 41 7 36]; description{5} = 'medium central spot electrodes';
%     electrodes{6} = [electrodes{5} 31 32 60 58 26 54 20 49 47 46 45 43 9 39 5 34]; description{6} = 'large central spot electrodes';
%     
%     elc_size_dist=[]; elc_ES_edges1=[];
%     for iElTopo=1:length(electrodes)
%         EEGtemp = pop_select( EEG, 'channel', electrodes{iElTopo});
%         
%         %find avalanches
%         elc_av_size_vec = [];
%         elc_sigma_vec = [];
%         for iEpoch = 1:EEG.trials
%             AvalancheResults = FindAvalanches(EEGtemp.data(:,:,iEpoch),EEGtemp.times,param.tau_vec(tau_optimal_inx),param.MaxMin,std_TH,1,0);
%             elc_av_size_vec = [elc_av_size_vec AvalancheResults.av_size_vec];
%             elc_sigma_vec = [elc_sigma_vec AvalancheResults.sigma_vec];
%         end 
%         elc_alpha = estimateParamML(param.ES.min,param.ES.max,'zeta',-param.optimal_alpha,elc_av_size_vec);
%         elc_sigma = mean(elc_sigma_vec); 
%         
%         % calculate avalanche distribution
%         [elc_size_dist(:,iElTopo), elc_ES_edges1(:,iElTopo)] = NormalizedDist(elc_av_size_vec,param.ES.edges); 
%         
%         figure;
%         subplot(1,2,1);topoplot([],EEGtemp.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEGtemp.chaninfo); title(description{iElTopo});
%         subplot(1,2,2); loglog(elc_ES_edges1(:,iElTopo),elc_size_dist(:,iElTopo),'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
%         set(gca,'XLim',[0 310]);xlabel('S');ylabel('P(S)');title(['\alpha = ' num2str(elc_alpha) ', \sigma = ' num2str(elc_sigma)]);
%     end
%      
%     figure; 
%     loglog(elc_ES_edges1,elc_size_dist,'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
%     set(gca,'XLim',[0 310]);xlabel('S');ylabel('P(S)');title(['\tau = ' num2str(1000*param.tau_vec(tau_optimal_inx)) ' ms']); legend(description);

end
