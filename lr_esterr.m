%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


function [press, estXs, betas] = lr_esterr(X,Y,metric)

n = size(X,1);
estXs = NaN(size(X,1),size(X,2));
betas = cell(size(X,1),size(X,2));

% predicted residual error
pre = zeros(1,n);


for j=1:n
    testdata = j;
    supervisors = setdiff(1:n, testdata);
    
    beta = lsqminnorm([ones(size(supervisors,2),1) Y(supervisors,:)], X(supervisors,:));
    betas(j,:) = mat2cell(beta,size(beta,1),ones(1,size(beta,2)));
    estXs(testdata,:) = [ones(size(testdata,2),1) Y(testdata,:)]*beta;
    
%     % another verion
%     for i=1:size(X,2)
%         beta = lsqminnorm([ones(size(supervisors,2),1) Y(supervisors,:)], X(supervisors,i));
%         betas{j,i} = beta;
%         estXs(testdata,i) = [ones(size(testdata,2),1) Y(testdata,:)]*beta;
%     end
    
    switch metric
        case 1
            pre(j) = sum((estXs(testdata,:)-X(testdata,:)).^2);
        case 2
            pre(j) =  1-dot(estXs(testdata,:),X(testdata,:))/(norm(estXs(testdata,:))*norm(X(testdata,:)));
        case 3
            corrr = 1-corrcoef(estXs(testdata,:)',X(testdata,:)');
            pre(j) = corrr(1,2);
        case 4
            pre(j) = norm(estXs(testdata,:)-X(testdata,:),1);
        case 5
            pre(j) = median(abs(estXs(testdata,:)-X(testdata,:))./(X(testdata,:)+1),2);
    end
    
end

% predicted residual error sum of squares
press = sum(pre);

end

