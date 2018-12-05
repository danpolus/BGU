function v = avs2vec(raster,grps,style)
%AVS2VEC represent avalanches as node#^2 vectors. The output is a count matrix 
%of prior/simultanious relations. the binary version of these vectors can be
%used as well. The edit distance is simply the absolute subtraction
%
%   INPUTS: binary raster (chan# x sample#)
%           a structure of electrode clicks (optional) which are lumped 
%           together using winner take all
%
%   OUTPUT: an array of avalanche topology vectors
if exist('grps','var')
    if ~exist('style','var')
       style = 0; 
    end
    raster = raster2grp(raster,grps,style);
end
avs = sum(raster);
[s,e] = enpoints2find(avs'>0);
if isempty(s)
    v = '';
    return
else
   v = nan(size(raster,1).^2,length(s));
   for n = 1 : length(s) 
       v(:,n) = av2vec(raster(:,s(n):e(n)));
   end
   v(:,sum(v)==0) = '';
end


function [s,e] = enpoints2find(msk)
%msk is assumed to be a column binary vector
dm = diff(msk);
e = find(dm==-1);
s = find(dm==1)+1;
if msk(1)
    s = [1;s];
end
if msk(end)
    e(end+1) = length(msk);
end

function v = av2vec(x)

[k,t] = size(x);
v = zeros(k);
alist = cellfun(@find,mat2cell(x,k,ones(1,t)),'uniformoutput',0);
for n = 1 : t
    inds = bsxfun(@plus,alist{n},(alist{n}'-1)*k);
    v(inds(:)) = v(inds(:)) + 1;
    if n < t
        inds = bsxfun(@plus,alist{n},(alist{n+1}'-1)*k);
        v(inds(:)) = v(inds(:)) + 1;
    end
end
v = v - diag(diag(v));
v = v(:);

function raster = raster2grp(raster,grps,style)

chanN = length(grps);
grpN = max(grps);
grpmat = zeros(grpN,chanN);
inds = sub2ind([grpN,chanN],grps,1:chanN);
grpmat(inds) = 1;
grpmat = bsxfun(@rdivide,grpmat,sum(grpmat,2));
raster = grpmat*raster;
if style == 1
    raster = raster > 0;
else
    raster = round(raster);
end


