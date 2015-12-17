%% Example with EEGLab STUDY
% EEGLAB is a popular tool for analyzing EEG data. Find out more about it
% here: http://sccn.ucsd.edu/wiki/EEGLAB
%
% The unbalancedPermutationGFP test is not well-integrated with EEGLAB, but
% with a little bit of coding, you can make it work.
%
% This example uses an EEGLAB STUDY structure to find the epoched EEG data
% on which to apply the unbalanced permutation test for GFP differences.
%
% The study used here is available from the EEGLAB wiki, and you can read
% about it here: http://sccn.ucsd.edu/wiki/Chapter_02:_STUDY_Creation
% The example data are approximately balanced, so the use of the
% unbalanced permutation test would not be necessary, but it still works as
% an example.
%
% Note EEGLAB has a lot of flexibility in how it actually stores the
% epoched data for a study. This script should work on the example data,
% but other studies that store their data differently will require
% potentially different treatment.
%
% Later on in the example, two corrections for multiple comparisons are
% also used: clusterCorrection and dmaxCorrection.

%% Setup
% -- nWorkers --
% unbalancedPermutationGFP will attempt to use the parallel toolbox.  If
% you do not have the parallel toolbox, this will be ignored. If you do
% have the parallel toolbox, but your system has limited memory available,
% set this to the number of subjects' data you can simultaneously load.
nWorkers = 12;
nPerm = 2000; % number of randomizations to use

study_dir = '~/matlab/gfp/example/5subjects/';
study_filename = 'n400clustedit.study';

%% Check that we can find EEGLab functions
if ~exist('pop_loadstudy', 'file')
    try
        eeglab('nogui');
    catch me
        e = error('UBpermGFP:NoEeglab',['I attempted to start EEGLAB, '...
            'but I failed. Make sure EEGLAB is either already running, '...
            'or is somewhere on MATLAB''s path.']);
        e.addCause(me);
        throw(e);
    end
end

%% Load the example study
% This loads all the data from the example study.  It uses about 550 MB in
% system memory, so it should be fine for systems with a few GB of memory
% free.
fprintf(1, 'Loading the study. This takes a while and throws warnings.\n');
starting_dir = pwd();
cd(study_dir);
[STUDY, ALLEEG] = pop_loadstudy('filepath', study_dir, ...
    'filename', study_filename);

%% Load the actual data for use with unbalancedPermutationGFP
FULLEEG = cell(size(STUDY.setind));
for iSet = 1:numel(FULLEEG),
    tmp = pop_loadset('eeg', ALLEEG(STUDY.setind(iSet)), ...
        'loadmode', 'all');
    tmp = pop_reref(tmp, []); % average reference the data!
    if iSet == 1,
        t = tmp.times;
    end
    FULLEEG{iSet} = tmp.data;
end
cd(starting_dir);
fprintf(1, 'Finished loading.\n');
%% Open a parallel pool
% 
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

%% Run the test
fprintf(1, 'Running the test.  This might take a while.\n');
tic
[p, AMean, BMean, sdists] = ...
    unbalancedPermutationGFP(FULLEEG(1,:), FULLEEG(2,:), nPerm);
fprintf(1, 'Done with permutations.\n');
toc

%% Make a figure showing the result
% Get middle 95% of permutation distribution
alpha = .05;
perm_d = AMean-BMean;
sorted_d = sort(perm_d, 2);

ci_lo = sorted_d(:, floor( (alpha/2)*nPerm));
ci_hi = sorted_d(:, ceil( (1-(alpha/2))*nPerm));

% some plotting variables
x = t;
patch_x = [x x(end:-1:1)];

% This figure shows the measured group mean GFP difference in the context
% of the middle 1-alpha of the permutation distribution.
figure('name', 'GFP Difference');
h1 = patch(patch_x, [ci_lo; ci_hi(end:-1:1)], [.6 .6 .6], ...
    'edgecolor', 'none');
hold on;
h2 = plot(x, perm_d(:, 1), 'b');
legend([h1 h2], {'Middle (1-\alpha) of permutation distribution', ...
    'Actual GFP difference'}, 'location', 'best');
xlabel('Time (ms)');
ylabel(sprintf('GFP difference (%s - %s)', STUDY.condition{:}));
title(STUDY.task);

%% Correct for multiple comparisons two fun and interesting ways
fprintf(1, 'Running corrections for multiple comparisons.\n');
tic
cluster_p = clusterCorrection(AMean, BMean, alpha);
dmax_p = dmaxCorrection(AMean, BMean);
fprintf(1, 'Done correcting for multiple comparisons.\n');
toc
% note: dmax_p returns p-values > 1.  These can safely be set to 1.
dmax_p(dmax_p>1) = 1;
%% Make a figure
figure('name', 'P-values');
h = plot(t, [p; cluster_p; dmax_p]);
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