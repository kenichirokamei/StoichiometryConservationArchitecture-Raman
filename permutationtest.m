%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


function [npress, significance, estXs, betas] = permutationtest(X,Y, npermutations, metric)

[npress, estXs, betas] = lr_esterr(X, Y, metric);

n = size(X,1);

if factorial(n)<npermutations
    npermutations = factorial(n)-1;
    rp = perms(1:n);
else
    rp = zeros(npermutations,n);
    for k=1:npermutations
        rp(k,:) = randperm(n);
    end
end

rpresses = zeros(1,npermutations);
for i=1:npermutations
    rpX = X(rp(i,:),:);
    switch metric
        case 1
            rpresses(i) = sum(sum((rpX-estXs).^2));
        case 2
            for j=1:n
                rpresses(i) = rpresses(i) + 1-dot(estXs(j,:),rpX(j,:))/(norm(estXs(j,:))*norm(rpX(j,:)));
            end
        case 3
            for j=1:n
                corrr = 1-corrcoef(estXs(j,:)',rpX(j,:)');
                rpresses(i) = rpresses(i) + corrr(1,2);
            end
        case 4
            for j=1:n
                rpresses(i) = rpresses(i) + norm(estXs(j,:)-rpX(j,:),1);
            end
        case 5
            rpresses(i) = sum(median(abs(estXs-rpX)./(rpX+1),2));
    end
end

% p-value correction according to 
% Phipson, B., and Smyth, G. K. (2010). Permutation p-values should never be zero: calculating exact
% p-values when permutations are randomly drawn. Stat. Appl. Genet. Molec. Biol. Volume 9, Issue
% 1, Article 39.

significance = (sum(rpresses<=npress)+1)/(npermutations+1);


end

