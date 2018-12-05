%
% calculate normalized similarity between 2 vectors of different length
% in case of 'jaccard' or 'correlation' similarity, slides short vector along the longer one
% return the maximal similarity found
%
function similarity = calc_sliding_similarity(v1, v2, bin_step_len, similarity_method)

similarity = [];

if strcmp(similarity_method,'levenshtein')
    similarity = calc_vectors_similarity(v1, v2, similarity_method);
    return;
end

%make sure v2 is longer than v1
if length(v1) > length(v2)
    v2_tmp = v2;
    v2 = v1;
    v1 = v2_tmp;
end

bins_v1 = length(v1)/bin_step_len;
bins_v2 = length(v2)/bin_step_len;
max_vec_len = max(length(v1),length(v2));

%       1 1 1
% 2 2 2 2 2
for i = 1:bins_v1-1
    v1_overlap = v1(1 : i*bin_step_len);
    v2_overlap = v2(end-i*bin_step_len+1 : end);
    similarity = [similarity calc_vectors_similarity(v1_overlap, v2_overlap, similarity_method)*length(v1_overlap)/max_vec_len];
end
%   1 1 1
% 2 2 2 2 2
for i = 1:bins_v2-bins_v1+1
    v1_overlap = v1;
    v2_overlap = v2(end-(bins_v1+i-1)*bin_step_len+1 : end-(i-1)*bin_step_len);
    similarity = [similarity calc_vectors_similarity(v1_overlap, v2_overlap, similarity_method)*length(v1_overlap)/max_vec_len];
end
% 1 1 1
%     2 2 2 2 2
for i = 1:bins_v1-1
    v1_overlap = v1(i*bin_step_len+1 : end);
    v2_overlap = v2(1 : (bins_v1-i)*bin_step_len);
    similarity = [similarity calc_vectors_similarity(v1_overlap, v2_overlap, similarity_method)*length(v1_overlap)/max_vec_len];
end

similarity = max(similarity);
