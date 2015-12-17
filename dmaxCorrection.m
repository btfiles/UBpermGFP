function dmax_p = dmaxCorrection(AMean, BMean)
% Correct for multiple comparisons using a max/min reference distribution.
% dmax_p = dmaxCorrection(AMean, BMean)
%
% INPUTS
% AMean and BMean are reference distributions of the group mean GFP for
% condition A and B, respectively. These are outputs 2 and 3 from
% unbalancedPermutationGFP(...).
%
% OUTPUT 
% dmax_p is a 1 by nsamp array with corrected p-values.  
%
% NOTES 
% This method builds a reference distribution of the maximal and minimal
% group mean GFP difference over the entire time course (rather than one
% reference distribution per time point, as with the uncorrected test).
% Then, the proportion of entries in the min distribution that are smaller
% than or equal to the actual value and the proportion of entries in the
% max distribution that are larger than or equal to the actual value are
% computed. The returned p-value is two times the smaller of these
% proportions.
%
% REFERENCE
% Blair, R. C., & Karniski, W. (1993). An alternative method for
% significance testing of waveform difference potentials. Psychophysiology,
% 30(5), 518?524.
%
% see also unbalancedPermutationGFP


mean_diff = AMean - BMean;
d_real = mean_diff(:, 1);
dmax_dist = max(mean_diff, [], 1);
dmin_dist = min(mean_diff, [], 1);
[ndmax, ndmin] = deal(zeros(size(d_real)));

for iT = 1:numel(d_real),
    ndmax(iT) = sum(dmax_dist >= d_real(iT));
    ndmin(iT) = sum(dmin_dist <= d_real(iT));
end

dmax_p = min([ndmax'; ndmin'])*2./size(AMean, 2);

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