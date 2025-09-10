%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


function plot_lda(U, dims, color, marker)

hold all

ngroups = size(U,1);

for i=1:ngroups
    plot(U{i}(:,dims(1)), U{i}(:,dims(2)),...
        'Marker',marker(mod(i-1,length(marker))+1),...
        'Color',color(mod(i-1,length(color))+1),...
        'LineStyle','none')
end

end

