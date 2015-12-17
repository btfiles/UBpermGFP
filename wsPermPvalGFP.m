function p = wsPermPvalGFP(distribution,probe)
% Compute p-values for permutation tests with null distributions that are
% not symmetrical about zero.
% p = wsPermPvalGFP(distribution,probe)
%
% INPUTS
% distribution is a nperm by n matrix
% probe is a 1 by n matrix
%
% OUTPUT
% p is the proportion of values in the distribution that are more extreme
% than the value in probe.  'More Extreme' is defined relative to the
% distribution. 
%
% NOTES
% This carries out a two-tailed test.  The minimum possible p-value to
% obtain is 2/nperm, so ensure that nperm is large enough to support the
% criterion p-value for your experiment.
%
%see also wsPermutationGFP unbalancedPermutationGFP

% validate probe
probe_size = size(probe);
if sum(probe_size ~= 1)~=1 && numel(probe) ~= 1,
    error('PermPGFP:BadProbe', ['Probe must have exactly 1 ', ...
        'non-singleton dimension (unless it is a scalar).']); 
end
probe = probe(:)'; % forces 1 by n shape.

% validate distribution
dist_size = size(distribution);
if dist_size(1) < 100,
    warning('PermPGFP:SmallDistribution', ['The reference ', ...
        'distribution is small (%d entries). This is fine for ', ...
        'debugging, but you should probably use a larger reference ', ...
        'distribution for hypothesis testing.'], dist_size(1));
end
assert(dist_size(2) == numel(probe), 'PermPGFP:SizeMismatch', ...
    ['The reference distribution and probe have incompatible sizes ', ...
     '(%d vs %d, respectively).'], dist_size(2), numel(probe));

% number of entries in the distribution that are greater than or equal
% to the probe
ngt = sum(bsxfun(@ge,distribution,probe)); 

% number of entries in the distribution that are less than or equal to the
% probe
nlt = sum(bsxfun(@le,distribution,probe));

% For a two-tailed test, compute number of more extreme entries
nme = 2*min([ngt; nlt]);

p = nme./dist_size(1);

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