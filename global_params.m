
function params_t = global_params()

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
params_t.raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well

params_t.similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well

params_t.compare_modes = {'Full', 'Len', 'ConcatLen'};%{'Len'};

params_t.test_set_percent = 0.1;

%minimal_similarity_threshold : to decide if matching cluster found and avoid anti-similarity
switch params_t.similarity_method
    case 'jaccard'
        params_t.minimal_similarity_threshold = 0.5;
    case 'correlation'
        params_t.minimal_similarity_threshold = 0;
    case 'levenshtein'
        params_t.minimal_similarity_threshold = 0.5;
end

nof_thresh_to_test = 10;
%params_t.max_cnt = 20; %for collapsing threshold. limit nof avalanches in descision
params_t.condition_descision_threshold = linspace(-1,4,nof_thresh_to_test); %static -1:4   step_av 1:5   step_stdv 1:5   contrast 0:1
params_t.condition_counter_limit = linspace(1,20,nof_thresh_to_test);
