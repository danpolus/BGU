%
%prepare artificially simulated vectors and Ids
%
%  inputs:
% MultiFileAchVecs - avalanche vectors 'prototype' 
% usedTauInfo - used taus info (optimal or others)
% saveFlg
%
%  outputs:
% SimMultiFileAchVecs - artificially simulated vectors and Ids
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SimMultiFileAchVecs = simulate_avalanche_vectors(MultiFileAchVecs, usedTauInfo, saveFlg)

max_Len = 3;
nof_ach_in_epoch = 5*max_Len;

%define avalanches of different lenght
sim_ach_vec{1} = zeros(1,60);
sim_ach_vec{1}(1:5) = 1;
sim_ach_vec{2} = zeros(2,60*2);
sim_ach_vec{2}(1,[6:10 71:75]) = 1;
sim_ach_vec{2}(2,[16:20 81:85]) = 1;
sim_ach_vec{3} = zeros(6,60*3);
sim_ach_vec{3}(1,[31:40 101:110 171:180]) = 1;
sim_ach_vec{3}(2,[31:40 111:120 161:170]) = 1;
sim_ach_vec{3}(3,[41:50 91:100 171:180]) = 1;
sim_ach_vec{3}(4,[41:50 111:120 151:160]) = 1;
sim_ach_vec{3}(5,[51:60 91:100 161:170]) = 1;
sim_ach_vec{3}(6,[51:60 101:110 151:160]) = 1;
%squareform(1 - pdist(sim_ach_vec{3},'jaccard'), 'tomatrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimMultiFileAchVecs = [];

iTau =  usedTauInfo.multi_files_tau_optimal_idx;

for iFile = 1:length(MultiFileAchVecs)
    SimMultiFileAchVecs{iFile}(iTau) = MultiFileAchVecs{iFile}(iTau);
    SimMultiFileAchVecs{iFile}(iTau).epochs_vecs = [];
    SimMultiFileAchVecs{iFile}(iTau).Id = [];
    SimMultiFileAchVecs{iFile}(iTau).Id{max_Len+1} = [];
    SimMultiFileAchVecs{iFile}(iTau).is_optimal_tau = true;   
    SimMultiFileAchVecs{iFile}(iTau).dataInfo.FileInfo.orig_fn = ['simulation_' SimMultiFileAchVecs{iFile}(iTau).dataInfo.FileInfo.orig_fn];
    
    for iEpoch = 1:length(MultiFileAchVecs{iFile}(iTau).epochs_vecs)
        for iAvalanche = 1:nof_ach_in_epoch
            switch mod(iAvalanche,max_Len)
                case 0
                    iLen = 3;
                    cond_idx = iFile;
                case 1
                    iLen = 1;
                    cond_idx = 1;
                case 2
                    iLen = 2;
                    cond_idx = str2num(MultiFileAchVecs{iFile}(iTau).file_id(13));
            end
            SimMultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch}(iAvalanche).vec = sim_ach_vec{iLen}(cond_idx,:);
            SimMultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch}(iAvalanche).first_bin_inx = iAvalanche;
            SimMultiFileAchVecs{iFile}(iTau).epochs_vecs{iEpoch}(iAvalanche).length_bins = iLen;
            
            SimMultiFileAchVecs{iFile}(iTau).Id{iLen} = [SimMultiFileAchVecs{iFile}(iTau).Id{iLen} {[MultiFileAchVecs{iFile}(iTau).file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')]}];
            SimMultiFileAchVecs{iFile}(iTau).Id{max_Len+1} = [SimMultiFileAchVecs{iFile}(iTau).Id{max_Len+1} {[MultiFileAchVecs{iFile}(iTau).file_id 'epc' num2str(iEpoch,'%03d') 'ach' num2str(iAvalanche,'%04d')]}];
        end       
    end
    
end

if saveFlg
    MultiFileAchVecs = SimMultiFileAchVecs;
    output_fp = [SimMultiFileAchVecs{1}(iTau).dataInfo.FileInfo.base_fp '2 avalanches\'];
    mkdir(output_fp);
    save([output_fp SimMultiFileAchVecs{1}(iTau).dataInfo.FileInfo.orig_fn '_avalanches.mat'],'MultiFileAchVecs','usedTauInfo');
end
