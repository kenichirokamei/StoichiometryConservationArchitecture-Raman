%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


%%%% Apply LDA to Raman spectra


%%% Randomly select 5/6th from data set

% 15 environmental conditions (BW25113)
conditions_used = 1:15; % alphabetical order

ngroups = length(conditions_used);

% For reproducibility
rng(1)

% Raman shift range
low = 2;
high = 600;

% Raman data smoothing
prefilter = @(x) savitzkyGolayFilt(x, 3, 0, 5, [], 2);
% Downloaded from 
% https://jp.mathworks.com/matlabcentral/fileexchange/30299-savitzky-golay-smooth-differentiation-filters-and-filter-application

% Raman data normalizing
postfilter = @(x) (x-mean(x,2))./std(x,0,2);

selraman = cell(ngroups,2);
remraman = cell(ngroups,2);    

selraman(:,1) = rawraman(conditions_used,1);
remraman(:,1) = rawraman(conditions_used,1);

for i=1:ngroups
    ramangroup = rawraman{conditions_used(i),2};
    groupsize = size(ramangroup,1);
    
    remn=floor(groupsize/6);
    seln = groupsize-remn;
    
    remindex = sort(randomsample_worep(seln,remn))'; % sampling uniformly at random, without replacement
    selindex = setdiff(1:groupsize, remindex);
    
    selraman{i,2} = postfilter(prefilter(rawraman{conditions_used(i),2}(selindex,low:high)));
    remraman{i,2} = postfilter(prefilter(rawraman{conditions_used(i),2}(remindex,low:high)));
end

% Spectra with cosmic ray spikes exist in 'selraman' created above. 
% The following code replaces the spectra with normal spectra in 'remraman'.
% Note that 'remraman' is not used in our analyses. 
    remramancell = 3;
    i=6; % Glucose42C
    j=1; % 1st cell
        selraman{i,2}(j,:) = remraman{i,2}(remramancell,:);
        remraman{i,2}(remramancell,:) = [];
    i=15; % Xylose
    j=23; % 23rd cell
        selraman{i,2}(j,:) = remraman{i,2}(remramancell,:);
        remraman{i,2}(remramancell,:) = [];


%%% Conduct LDA

U = lda(selraman(:,2), selraman(:,1), 98);

% calculate mean of each group
meanU = NaN(ngroups, ngroups-1);
for i=1:ngroups
    meanU(i,:) = mean(U{i,2}, 1);
end





