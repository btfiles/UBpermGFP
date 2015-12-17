%% UBpermGFP Example: file list syntax
% example.m
%
% This file goes through an example of doing an unbalanced permutation test
% for differences in group mean global field power. This uses the list
% syntax for unbalancedPermutationGFP(...), which should work for systems
% with less memory. (See the note about nWorkers, below)
%
% This produces a figure summarizing the GFP difference and a second figure
% showing p-values as corrected for multiple comparisons using the
% cluster-based and dmax-based corrections.

%% Setup
global VERBOSE
VERBOSE = true;

% -- nWorkers --
% unbalancedPermutationGFP will attempt to use the parallel toolbox.  If
% you do not have the parallel toolbox, this will be ignored. If you do
% have the parallel toolbox, but your system has limited memory available,
% set this to the number of subjects' data you can simultaneously load.
nWorkers = 12;

nPerm = 2000; % number of randomizations to use

field = 'd'; % the name of the variable inside your .mat file that has 
             % epoched data
alpha = .05; % criterion p-value (used for generating a figure)

%% Create file lists
% This simply creates two cell arrays, A_list and B_list, with each cell
% containing a fully-qualified filename.  These files/paths probably do not
% exist on your system, so you will need to edit this section to get the
% code to run.
epoch_dir = '~/example/epochs/';
epoch_pattern_A = 'S%02d_TargetEpochs.mat';
epoch_pattern_B = 'S%02d_BackgroundEpochs.mat';
sid_list = 1:13;

[A_list, B_list] = deal(cell(1, numel(sid_list)));
for iSub = 1:numel(sid_list),
    sid = sid_list(iSub);
    A_list{iSub} = fullfile(epoch_dir, sprintf(epoch_pattern_A, sid));
    B_list{iSub} = fullfile(epoch_dir, sprintf(epoch_pattern_B, sid));
end

%% Open a parallel pool
% Should work without parallel computing toolbox
if license('test', 'Distrib_Computing_Toolbox') && exist('gcp', 'file'),
    if isempty(gcp('nocreate')),
        parpool(nWorkers);
    end
elseif nWorkers > 1,
    warning('UBpermGFP:NonParallel', ['The parallel computing toolbox ',...
        'does not appear to be available, but you requested %d ',...
        'workers. Falling back to 1 worker.'], nWorkers);
    nWorkers = 1;
end
%% Run the file list version
% outputs p and sdists are not used in this example, but are shown for
% syntax illustration.
tic
fprintf(1, 'Running test...\n');
[p, AMean, BMean, sdists] = ...
    unbalancedPermutationGFP(A_list, B_list, nPerm, field);
fprintf(1, 'Done running test.\n');
toc
%% Make a figure showing the result
% Get middle 95% of permutation distribution
perm_d = AMean-BMean;
sorted_d = sort(perm_d, 2);

ci_lo = sorted_d(:, floor( (alpha/2)*nPerm));
ci_hi = sorted_d(:, ceil( (1-(alpha/2))*nPerm));

% some plotting variables
x = 1:size(perm_d,1);
patch_x = [x x(end:-1:1)];

% This figure shows the measured group mean GFP difference in the context
% of the middle 1-alpha of the permutation distribution.
figure;
h1 = patch(patch_x, [ci_lo; ci_hi(end:-1:1)], [.6 .6 .6], ...
    'edgecolor', 'none');
hold on;
h2 = plot(x, perm_d(:, 1), 'b');
legend([h1 h2], {'Middle (1-\alpha) of permutation distribution', ...
    'Actual GFP difference'});
xlabel('Time (samples)');
ylabel('GFP difference (A-B)');

%% Correct for multiple comparisons two fun and interesting ways
fprintf(1, 'Running corrections for multiple comparisons...\n');
tic
cluster_p = clusterCorrection(AMean, BMean, alpha);
dmax_p = dmaxCorrection(AMean, BMean);
fprintf(1, 'Done with multiple comparisons correction.\n');
toc
% note: dmax_p returns p-values > 1.  These can safely be set to 1.
dmax_p(dmax_p>1) = 1;
%% Make a figure
figure('name', 'P-values');
h = plot(x, [p; cluster_p; dmax_p]);
ylim([0 .15]);
% set(gca,'YScale', 'log');
hold on;
plot(xlim(), alpha.*[1 1],'k--');
legend(h, {'Uncorrected', 'Cluster correction', 'd-max correction'}, ...
    'location', 'best');
xlabel('time (ms)');
ylabel('p-value');

% Change Log
%   2015-12-04 Initial Release - BTF
%
%
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