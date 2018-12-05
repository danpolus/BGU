%
% calculate normalized similarity between 2 vectors
%
% vectors have same length for jaccard or correlation (not neccessary for levenshtein)
% levenshtein distance with MAX (not SUM) non zero normalization is used
%
function similarity = calc_vectors_similarity(v1, v2, similarity_method)

similarity = [];

switch similarity_method
    case 'jaccard'
        %jaccard (boolean) distance:
        %pdist range 0-:-1
        %good for binary vectors of same length
        %works for double too but don't consider values differences, just the fact they're different
        similarity = 1 - pdist([v1; v2], 'jaccard');
    case 'correlation'
        %correlation distance:
        %pdist range 0-:-2 (supports anti-correlations)
        %good for double vectors of same length
        non_zero_idx = v1~=0 | v2~=0; %exclude mutual zeroes
        if std(v1(non_zero_idx)) == 0 || std(v2(non_zero_idx)) == 0
            similarity = 1 - pdist([v1(non_zero_idx) 0; v2(non_zero_idx) 0], 'correlation'); %add 0 to vectors to avoid std=0
        else
            similarity = 1 - pdist([v1(non_zero_idx); v2(non_zero_idx)], 'correlation');
        end
    case 'levenshtein'
        %levenshtein distance:
        %lev range 0-:-max(len(v1),len(v2)).  range after normalization typicaly 0-:-1, but can be >1 ("anti-correlation")
        %good for binary vectors
        %works for double too but don't consider values differences, just the fact they're different
        %levenshtein distance:
        lev_dist = lev(v1,v2);
        if length(v1)==length(v2)
            %max non zero normalization - less forgivable. better when copmaring avalanches
            similarity = 1 - lev_dist / max(sum(v1~=0), sum(v2~=0));
            %%sum non zero normalization - more forgivable
            %similarity = 1 - lev_dist / sum([v1 v2]~= 0);
        else
            lev_dist = lev_dist - abs(length(v1)-length(v2));
            %max non zero normalization - less forgivable. better when copmaring avalanches
            lev_dist = lev_dist / max(sum(v1~=0), sum(v2~=0));
            %%sum non zero normalization - more forgivable
            %lev_dist = lev_dist / sum([v1 v2]~=0);
            similarity = 1 - lev_dist * max(length(v1),length(v2))/min(length(v1),length(v2));  
        end
end
