%
%Step 5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function predict_conditions()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select clustering results files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
    load([fp files{iFile}],'ClusteringData');
    load([fp [files{iFile}(1:end-13) '.mat']],'MultiFileAchVecs','SimilarityMat','TestingSet');
    
    PredictionResults = [];
    for iTau = 1:length(ClusteringData)
        
        if isempty(ClusteringData(iTau).tau)
            continue;
        end
        
        for iCond=1:length(TestingSet(iTau).EpochIds)
            %v1 = MultiFileAchVecs{str2num(SimilarityMat(iTau).Id{m}(21:23))}(iTau).epochs_vecs{str2num(SimilarityMat(iTau).Id{m}(27:29))}(str2num(SimilarityMat(iTau).Id{m}(33:36))).vec;
            clst=find_vector_cluster(v);
            [cond, certainty] = bayessian_accumulative_predictor(clst);
        end
        
        
    end
    
    save([fp files{iFile}(1:end-4) '_prediction.mat'],'PredictionResults');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cond, certainty] =  bayessian_accumulative_predictor(clst)

threshold_ = 0.5;
%%%%%%%%%%%%%%%%%

persistent accumulator;
%exponential primacy weighting?
%static thresholds? distribution depandent?
%dynamic thresholds? accumulator dynamic range and std depandent?
%threshold collapses over time? as function of nof samples?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_num = find_vector_cluster(v)
%for each cluster calculate average similarity between v and cluster's vectors
%return the cluster with maximal similarity

minimal_distance_threshold = 0.7;
%%%%%%%%%%%%%%%%

