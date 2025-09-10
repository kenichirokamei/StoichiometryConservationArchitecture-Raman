%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


function [U, nPCs, X] = lda(data, labels, percent)

ngroups = size(data,1);

X = [];

groupnames = {};

for i=1:ngroups
    groupsize = size(data{i},1);
    X(end+1:end+groupsize,:) = data{i};
    groupnames(end+1:end+groupsize)=labels(i);
end

% PCA
Xcentered = X-mean(X,1);
[~,Xcentered_S,Xcentered_V] = svd(Xcentered,"econ");
    % For reproducibility, make the sign of the element with the largest absolute value in each
    % column positive.
    % (Signs of singular vectors that 'svd' function of Matlab outputs can be different 
    % across different machines and releases of Matlab.)
    [~,rowind_eachcol] = max(abs(Xcentered_V),[],1);
    Xcentered_V = Xcentered_V .* sign(Xcentered_V(sub2ind(size(Xcentered_V),rowind_eachcol,1:size(Xcentered_V,2))));
ramanscore = Xcentered*Xcentered_V;
ramanexplained = diag(Xcentered_S).^2/sum(diag(Xcentered_S).^2)*100;

% use PCs that in total explain 'percent'% of variation
hoge = find(cumsum(ramanexplained)>percent);
nPCs = hoge(1);

% LDA
C_I = zeros(nPCs,nPCs); % cov within groups
M_I = NaN(ngroups,nPCs); % mean within groups
for i=1:ngroups
    index = strcmp(labels{i}, groupnames);
    M_I(i,:) = mean(ramanscore(index,1:nPCs),1);
    C_I = C_I + cov(ramanscore(index,1:nPCs));
end
C_I = C_I/ngroups;
C_B = cov(M_I); % cov between groups
[Vmat,eee] = eigs(C_B, C_I, min(nPCs, ngroups-1));

[~, indexeee] = sort(diag(eee));
Vmat = Vmat(:,indexeee(end:-1:1));
V = normalize(Vmat,1,'norm');

Utmp = ramanscore(:,1:nPCs)*V;

U = cell(ngroups,2);
for i=1:ngroups
    index = strcmp(labels{i}, groupnames);
    U{i,2}=Utmp(index,:);    
end

U(:,1) = labels;

end

