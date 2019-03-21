
function params_t = global_params()

%the raster matrix is not necessary binary. it counts events per bin, so it's values can be >1. it can also contain amplitudes
%for now, make it binary (OTHER OPTIONS NEVER TESTED OR DEBUGED)
params_t.raster_input_type = 'binary'; % 'binary' 'events_cnt' 'amplitude' -> set to appear in folder name as well

params_t.similarity_method = 'jaccard'; % 'jaccard' 'correlation' 'levenshtein'  -> set to appear in folder name as well

params_t.compare_modes = {'Full', 'Len', 'ConcatLen'};%{'Len'};

params_t.test_set_percent = 0.1;

%minimal_similarity_threshold : to decide if matching cluster found. depends on similarity_method
switch params_t.similarity_method
    case 'jaccard'
        params_t.minimal_similarity_threshold = 0.5; %below 0.5 it's anti-similarity
    case 'correlation'
        params_t.minimal_similarity_threshold = 0.5;
    case 'levenshtein'
        params_t.minimal_similarity_threshold = 0.5;
end
params_t.condition_descision_threshold = 5; %threshold for accumulative predictor to decide
params_t.accum_sample_limit = 10; %max nof avalanches in accumulator for SampLimitAccum mode
