%
%Step 2 : extract_avalanches - extract avalanches from EEGlab sets
%
%  inputs:
% EEGSets - cleaned data separated into words and conditions, as EEGlab sets
% plotFlg
%
%  outputs:
% AvalancheFileDataSets - avalanches extracted from each EEGlab set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AvalancheFileDataSets = extract_avalanches(EEGSets, plotFlg)

params_t = global_params();

%prepare global params for FindAvalanches
global param
param.Fs = EEGSets(1).srate;
%Event Size - range of avalanche sizes
param.ES.min = 1;
param.ES.max = (floor(EEGSets(1).nbchan*2/100))*100;
param.ES.edges = unique(ceil(logspace(0,log10(param.ES.max),25))); % log spaced bins
%times for raster plot
param.t1 = 4; % sec
param.t2 = 8; % sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AvalancheFileDataSets = [];

subj_scenario_bad_channels = [];

%calculate mean and std of the whole scenario for each channel
if strcmp(params_t.zscore_mode, 'scenario')
    EEG_data_reshaped = [];
    for iEegSets = 1:length(EEGSets)
        EEG = EEGSets(iEegSets);
        bad_epoch_chan = boolean(EEG.reject_hstr.rejglobalE(:,~EEG.reject_hstr.rejglobal));
        for iEpoch=1:EEG.trials
            EEG.data(bad_epoch_chan(:,iEpoch),:,iEpoch) = NaN; %set bad epoch channels to NaN
        end
        EEG_data_reshaped = [EEG_data_reshaped reshape(EEG.data,size(EEG.data,1),[])];
        subj_scenario_bad_channels = [subj_scenario_bad_channels EEG.bad_channels];
    end
    subj_mean = mean(EEG_data_reshaped,2,'omitnan');
    subj_std = std(EEG_data_reshaped,[],2,'omitnan');
    subj_scenario_bad_channels = unique(subj_scenario_bad_channels);
    clear('EEG_data_reshaped');
end

%extract avalanches
for iEegSets = 1:length(EEGSets)
    EEG = EEGSets(iEegSets);
    
    %save data specific params
    dataInfo.fs = EEG.srate;
    dataInfo.tau_vec = 1/EEG.srate:1/EEG.srate:params_t.max_tau_sec; % array with delta t values for avalanche analysis
    dataInfo.FileInfo = EEG.FileInfo;
    dataInfo.ES = param.ES;
    
    %handle bad channels and epochs
    bad_epoch_chan = boolean(EEG.reject_hstr.rejglobalE(:,~EEG.reject_hstr.rejglobal));
    if params_t.zero_subj_scenario_bad_channels_flg
        bad_epoch_chan(subj_scenario_bad_channels,:) = true;
    end
    for iEpoch=1:EEG.trials
        EEG.data(bad_epoch_chan(:,iEpoch),:,iEpoch) = NaN; %set bad epoch channels to NaN
    end
    
    %zscore
    if strcmp(params_t.zscore_mode, 'file')
        subj_mean = mean(reshape(EEG.data,size(EEG.data,1),[]),2,'omitnan');
        subj_std = std(reshape(EEG.data,size(EEG.data,1),[]),[],2,'omitnan');
    end
    if strcmp(params_t.zscore_mode, 'scenario') || strcmp(params_t.zscore_mode, 'file')
        %EEG.data = zscore(EEG.data,0,[2 3]); %for file: should work in matlab 2019
        subj_mean_mat = repmat(subj_mean,[1,size(EEG.data,2),size(EEG.data,3)]);
        subj_std_mat = repmat(subj_std,[1,size(EEG.data,2),size(EEG.data,3)]);
        EEG.data = (EEG.data - subj_mean_mat)./ subj_std_mat;
    end
    if strcmp(params_t.zscore_mode, 'epoch')
        EEG.data = normalize(EEG.data,2,'zscore');
    end
    EEG.data(isinf(EEG.data) | isnan(EEG.data)) = 0; %set bad epoch channels to 0
    
    %plot
    if plotFlg
        EEG_plot = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -params_t.std_TH, params_t.std_TH, EEG.xmin, EEG.xmax, 0, 0);
        EEG_plot = eeg_rejsuperpose(EEG_plot, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG_plot.reject.rejmanual = EEG_plot.reject.rejglobal;
        EEG_plot.reject.rejmanualE = EEG_plot.reject.rejglobalE;
        pop_eegplot(EEG_plot, 1, 0, 0, [], 'srate',EEG_plot.srate, 'winlength',10, 'spacing', 5, 'eloc_file', []);
    end
    
    %initialize results vector
    all_epochs = [];
    for iTau = 1:length(dataInfo.tau_vec)
        all_epochs(iTau).av_raster_epochs = [];
        all_epochs(iTau).av_raster = [];
        all_epochs(iTau).av_size_vec = [];
        all_epochs(iTau).av_dur_vec = [];
        all_epochs(iTau).sigma_vec = [];
        all_epochs(iTau).tau = dataInfo.tau_vec(iTau);
    end
    
    %update params
    param.Fs = EEG.srate;
    param.ES.max = (floor(EEG.nbchan*2/100))*100;
    param.ES.edges = unique(ceil(logspace(0,log10(param.ES.max),25)));
    dataInfo.ES = param.ES;
    
    %find avalanches
    for iEpoch = 1:EEG.trials
        AvalancheResults = FindAvalanches(EEG.data(:,:,iEpoch),EEG.times,dataInfo.tau_vec,'maxmin',params_t.std_TH,1,0);
        for iTau = 1:length(dataInfo.tau_vec)
            all_epochs(iTau).av_raster_epochs(:,:,iEpoch) = AvalancheResults(iTau).av_raster;
            all_epochs(iTau).av_raster = [all_epochs(iTau).av_raster AvalancheResults(iTau).av_raster];
            all_epochs(iTau).av_size_vec = [all_epochs(iTau).av_size_vec AvalancheResults(iTau).av_size_vec];
            all_epochs(iTau).av_dur_vec = [all_epochs(iTau).av_dur_vec AvalancheResults(iTau).av_dur_vec];
            all_epochs(iTau).sigma_vec = [all_epochs(iTau).sigma_vec AvalancheResults(iTau).sigma_vec];
        end
    end
    
    for iTau = 1:length(dataInfo.tau_vec)
        all_epochs(iTau).alpha = estimateParamML(param.ES.min,param.ES.max,'zeta',-params_t.optimal_alpha,all_epochs(iTau).av_size_vec);
        all_epochs(iTau).sigma = mean(all_epochs(iTau).sigma_vec);
        
        if plotFlg % rate of events in each channel
            R(:,iTau) = sum(all_epochs(iTau).av_raster,2)/(size(all_epochs(iTau).av_raster,2)/EEG.srate);
        end
    end
    
    AvalancheFileDataSets(iEegSets).all_epochs = all_epochs;
    AvalancheFileDataSets(iEegSets).dataInfo = dataInfo;
    
    %plot
    if plotFlg
        [~,tau_optimal_inx] = min(sqrt(([all_epochs.sigma] - params_t.optimal_sigma).^2 + (-[all_epochs.alpha] - params_t.optimal_alpha).^2));
        
        xx = 1:300; yy = xx.^(params_t.optimal_alpha);
        xSize = 7; ySize = 9.5; xLeft = (8.5 - xSize)/2; yTop = (11 - ySize)/2;%Graphic parameters
        for iTau = 1:length(dataInfo.tau_vec)
            
            if iTau~=tau_optimal_inx
                continue;
            end
            
            figure('Name',[EEG.setname '  \tau=' num2str(iTau) 'dt']);
            set(gcf,'Color','w');
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
            set(gcf,'Position',[570 40 xSize*75 ySize*75])
            
            subplot(2,3,1);histogram(all_epochs(iTau).av_size_vec,max(all_epochs(iTau).av_size_vec));title(['Avalanche Size  \tau = ' num2str(dataInfo.tau_vec(iTau))]);
            subplot(2,3,4);histogram(all_epochs(iTau).av_dur_vec,max(all_epochs(iTau).av_dur_vec));title(['Avalanche Duration  \tau = ' num2str(dataInfo.tau_vec(iTau))]);
            
            subplot(2,3,2); RasterPlot(all_epochs(iTau).av_raster,param.Fs,dataInfo.tau_vec(iTau));title(['\tau = ' num2str(dataInfo.tau_vec(iTau)) ' - Raster']); drawnow
            subplot(2,3,3); RasterPlot(all_epochs(iTau).av_raster,param.Fs,dataInfo.tau_vec(iTau),[param.t1 param.t2]);title('zoom in'); drawnow
            subplot(2,3,5); stem(R(:,iTau));xlabel('Channel #');ylabel('Events/sec');title('Event rate in each channel');set(gca,'XLim',[1 size(all_epochs(iTau).av_raster,1)]);
            
            % calculate avalanche distribution
            [size_dist, ES_edges1] = NormalizedDist(all_epochs(iTau).av_size_vec,param.ES.edges);
            subplot(2,3,6); loglog(ES_edges1,size_dist,'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
            set(gca,'XLim',[0 310]);
            xlabel('S');ylabel('P(S)');
            title(['\tau = ' num2str(1000*dataInfo.tau_vec(iTau)) ' ms' ', \alpha = ' num2str(all_epochs(iTau).alpha) ', \sigma = ' num2str(all_epochs(iTau).sigma)]);
            
        end
        
        figure('Name',EEG.setname);scatter([all_epochs.sigma],-[all_epochs.alpha],'x');xlabel('\sigma');ylabel('\alpha');title('\alpha = F(\sigma)  labels: \tau values [sec]');
        hold on
        scatter([all_epochs(tau_optimal_inx).sigma],-[all_epochs(tau_optimal_inx).alpha],'dk','filled');
        text([all_epochs.sigma]+0.01,-[all_epochs.alpha],strsplit(num2str(dataInfo.tau_vec)));
        plot([min(params_t.optimal_sigma,min([all_epochs.sigma]))-0.1, max(params_t.optimal_sigma,max([all_epochs.sigma]))+0.1],[params_t.optimal_alpha,params_t.optimal_alpha],'--m',...
            [params_t.optimal_sigma,params_t.optimal_sigma],[min(params_t.optimal_alpha,min(-[all_epochs.alpha]))-0.1, max(params_t.optimal_alpha,max(-[all_epochs.alpha]))+0.1],'--m');
        hold off
    end
    
    
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
%             AvalancheResults = FindAvalanches(EEGtemp.data(:,:,iEpoch),EEGtemp.times,dataInfo.tau_vec(tau_optimal_inx),'maxmin',params_t.std_TH,1,0);
%             elc_av_size_vec = [elc_av_size_vec AvalancheResults.av_size_vec];
%             elc_sigma_vec = [elc_sigma_vec AvalancheResults.sigma_vec];
%         end
%         elc_alpha = estimateParamML(param.ES.min,param.ES.max,'zeta',-params_t.optimal_alpha,elc_av_size_vec);
%         elc_sigma = mean(elc_sigma_vec);
%         
%         % calculate avalanche distribution
%         [elc_size_dist(:,iElTopo), elc_ES_edges1(:,iElTopo)] = NormalizedDist(elc_av_size_vec,param.ES.edges);
%         
%         if plotFlg
%             figure;
%             subplot(1,2,1);topoplot([],EEGtemp.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEGtemp.chaninfo); title(description{iElTopo});
%             subplot(1,2,2); loglog(elc_ES_edges1(:,iElTopo),elc_size_dist(:,iElTopo),'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
%             set(gca,'XLim',[0 310]);xlabel('S');ylabel('P(S)');title(['\alpha = ' num2str(elc_alpha) ', \sigma = ' num2str(elc_sigma)]);
%         end
%     end
%     
%     if plotFlg
%         figure;
%         loglog(elc_ES_edges1,elc_size_dist,'LineWidth',2); hold on; loglog(xx,yy,'k--','LineWidth',3); hold off;
%         set(gca,'XLim',[0 310]);xlabel('S');ylabel('P(S)');title(['\tau = ' num2str(1000*dataInfo.tau_vec(tau_optimal_inx)) ' ms']); legend(description);
%     end
    
end
