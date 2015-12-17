function [pval,aMeanDist,bMeanDist, sDists] = ...
    unbalancedPermutationGFP(A, B, nPerm, flds)
%UNBALANCEDPERMUTATIONGFP Runs the unbalanced permutation test. 
% [pval,AMean, BMean, sDists] = unbalancedPermutationGFP(A,B,nPerm)
%
% INPUTS
%
% A and B are 1xN cell arrays of data matrices.  A{n} is the collection of
% trials/epochs labeled A for subject n, B{n} is the collection of
% trials/epochs labeled B for subject n. Data matrices should be nbchan by
% nsamp by ntrial (following the EEGLab convention for epoched data). It is
% critical that the data have already been average referenced.
%
% nPerm is the number of permutations/randomizations requested.  This
% should be a pretty big number; under fairly modest assumptions the
% ceiling on this number is too big to consider, so the upper limit is
% really about how long you want to wait for an answer balanced against the
% p-value resolution you're going to want.  The code complains if this is
% less than 100, but something like 2000 is probably okay for a lot of
% purposes.
%
% OUTPUTS 
%
% pval is a two-tailed p-value for obtaining a mean GFP difference as
% extreme as or more extreme than the one actually obtained given the null
% hypothesis that the labels for data (A or B) are random or arbitrary.
%
% aMeanDist, bMeanDist are the group mean GFPs in the null distribution.
% These are nPerm x nSamp, where nSamp is the number of EEG samples per
% epoch.
%
% sDists is a 2 by n cell array with each cell containing the
% single-subject permutation distribution for each subject. Condition A' is
% in row 1, condition B' is in row 2.
%
% NOTE
% The data-as-labeled are always added to the permutation distribution as
% the first 'permutation'.
%
% LOW MEMORY SYSTEMS or really, SYSTEMS WITHOUT GOBS OF FREE MEMORY
%
% If you do not have enough system memory to keep all the data for your
% experiment loaded simultaneously, an alternative syntax is provided. This
% should not even be slower than the first syntax, but might be less
% convenient. Note, you still need enough memory to load data for both
% conditions from a sinlge subject.
%
% [...] = unbalancedPermutationGFP(A, B, nPerm, fields)
% where A and B are cell arrays of strings containing file names. The file
% names must either be fully qualified in a way the system understands or
% point to files already on the matlab path. As above, they are both 1xN
% for N subjects.  Argument fields is either a string or cell array of
% strings for specifying how to get the data out of your file.
% L = load(A{idx}, '-mat'); 
% data_a_idx = L.(fields); 
% OR
% data_a_idx = L.(fields{1}).(fields{2}).(fields{n})
%
% For example, if your epoched data are in .mat files in a variable called
% 'data', call [...] = unbalancedPermutationGFP(A, B, 2000, 'data').
% If your data are in EEGLAB .set files (and you are not using the option
% to store data and headers in separate files), then the call would be
% [...] = unbalancedPermutationGFP(A, B, 2000, {'EEG', 'data'}).
%
% This does not cover every use case for data storage, but this function
% should be pretty straightforward to modify to cover your storage case.
% Or, if you were feeling really ambitious, you could overload matlab's
% inbuilt load function.
%
%see also wsPermutationGFP wsPermPvalGFP

global VERBOSE
if ~exist('VERBOSE', 'var')
    verbose = false;
else
    verbose = VERBOSE;
end

% Validate inputs
assert(numel(A)==numel(B), ...
    'Inputs A and B must have the same number of elements.');

if nPerm < 100,
    warning('UBGFP:LowPermutations', ['You have requested a relatively',...
        ' small number of permutations (%d). This is fine for ', ...
        'debugging, but for hypothesis testing you should use more.'], ...
        nPerm);
end

if ~exist('flds', 'var'),
    flds = [];
end
lomem = ischar(A{1});
% Validate input files
if lomem,
    if verbose,
        disp('Using low-memory syntax for unbalanced permutation test.');
    end
    for iA = 1:numel(A),
        assert(exist(A{iA}, 'file')>1,'UBGFP:LowMem:FileNotFound', ...
            'File %s could not be found.', A{iA});
    end
    for iB = 1:numel(B),
        assert(exist(B{iB}, 'file')>1,'UBGFP:LowMem:FileNotFound', ...
            'File %s could not be found.', B{iB});
    end
    % Wrap scalar fields in a cell
    if ~iscell(flds),
        flds = {flds};
    end
end

% Now get to work
sDists = cell(2,length(A));
parfor iSub = 1:length(A),
    if verbose,
        display(['Working on sub ' num2str(iSub)]);
    end
    if lomem,
        dA = load_data(A{iSub}, flds);
        dB = load_data(B{iSub}, flds);
        [aPerm, bPerm] = wsPermutationGFP(dA, dB, nPerm);
    else
        [aPerm, bPerm] = wsPermutationGFP(A{iSub}, B{iSub}, nPerm);
    end
    sDists(:,iSub) = {squeeze(aPerm); squeeze(bPerm)};
end
if verbose,
    display('All workers done.');
end

aMeanDist = mean(cat(3,sDists{1,:}),3);
bMeanDist = mean(cat(3,sDists{2,:}),3);

dMeanDist = aMeanDist-bMeanDist;
pval = wsPermPvalGFP(dMeanDist',dMeanDist(:,1)');

end % end unbalancedPermutationGFP

% Helper function
function d = load_data(fn, fields)
d = load(fn, '-mat');
% dig through possibly nested structs to get to the juicy data at the end
for idx = 1:numel(fields),
    d = d.(fields{idx});
end
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