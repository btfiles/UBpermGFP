function clust_p = clusterCorrection(AMean, BMean, alpha)
%Correct for multiple comparisons using the cluster size correction.
%clust_p = clusterCorrection(AMean, BMean, alpha)
%
%AMean and BMean are reference distributions of the group mean GFP for
%condition A and B, respectively. These are outputs 2 and 3 from
%unbalancedPermutationGFP(...).
%
%alpha is the (uncorrected) criterion p-value for inclusion in a cluster.
%
%Reference: Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S.,
%Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, by
%theory and permutation, for a difference between two groups of structural
%MR images of the brain. IEEE Trans Med Imaging, 18(1), 32?42.
%
%Note: This function uses a parfor loop. The mathworks say it should work
%(albeit not in parallel) even without the parallel toolbox.
%
%see also wsPermPvalGFP unbalancedPermutationGFP

% Cluster-based
md = AMean' - BMean';
size_mx = zeros(size(md,1), 1);
[real_clust_size, real_p, real_onsets, real_offsets] = ...
    deal([]); % trick parfor into keeping these

% Check that we have at least some clusters:
p = wsPermPvalGFP(md, md(1, :));
if ~any(p<alpha),
    % If not, we are done!
    clust_p = ones(size(p));
    return;
end

% Iterate over each randomization/permutation to build a null distribution
% of cluster sizes
parfor i = 1:size(md, 1),
    
    % Compute p-values for this row relative to all other rows
    p = wsPermPvalGFP(md, md(i, :));
    s = p < alpha;
    if any(s)
        % find clusters of consecutive p-values < alpha.
        [clust_onset, clust_offset] = deal(false(size(s)));
        clust_onset(2:end) = s(2:end) & ~s(1:(end-1));
        clust_onset(1) = s(1);
        
        clust_offset(1:(end-1)) = s(1:(end-1)) & ~ s(2:end);
        clust_offset(end) = s(end);
        
        onset_idx = find(clust_onset);
        offset_idx = find(clust_offset);
        
        % cluster size and max for this row
        clust_size = offset_idx-onset_idx+1;
        size_mx(i) = max(clust_size);
        
        % hang onto values for the first row (the real data)
        if i==1,
            real_clust_size(i,:) = clust_size;
            real_p(i,:) = p;
            real_onsets(i,:) = onset_idx;
            real_offsets(i,:) = offset_idx;
        end
    end
end

clust_p = ones(size(real_p));

for iReal = 1:numel(real_clust_size),
    p = mean(size_mx>=real_clust_size(iReal));
    clust_p(real_onsets(iReal):real_offsets(iReal)) = p;
end

% Change Log
%   2015-12-04 Initial Release - BTF

% Copyright notice
%    Copyright 2015 Benjamin T. Files
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
%    implied. See the License for the specific language governing
%    permissions and limitations under the License.