

function y = randomsample_worep(n,k)

    y = [];
    nvals = 0;
    while nvals<k
        y(end+1:end+k-nvals) = randi(n,1,k-nvals);
        y = unique(y);
        nvals = length(y);
    end
    y = y(randperm(k));

end

