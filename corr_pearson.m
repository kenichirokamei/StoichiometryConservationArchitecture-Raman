%%% Copyright 2023 Ken-ichiro F. Kamei %%%


function R = corr_pearson(X,form)
    
    ncol = size(X,2);
    R = NaN(ncol,ncol);
    comb = nchoosek(1:ncol,2);
    for i=1:size(comb,1)
        tmp = corrcoef(X(:,comb(i,:)));
        R(comb(i,1),comb(i,2)) = tmp(1,2);
    end
    
    switch form
        case 'matrix'
            R = triu(R,1)+triu(R,1)'+eye(ncol);
        case 'vector'
            R = R(logical(triu(ones(ncol,ncol),1)));
        otherwise
            error("Invalid parameter name: "+string(form))
    end
    
end

