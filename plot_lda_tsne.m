%%% Copyright 2018-2023 Ken-ichiro F. Kamei %%%


function plot_lda_tsne(U, color, marker)

U_fortsne = [];
U_names_fortsne = [];
for i=1:size(U,1)
    U_fortsne = [U_fortsne; U{i,2}];
    U_names_fortsne = [U_names_fortsne; repmat(string(U{i,1}), size(U{i,2},1), 1)];
end

tsneresult = tsne(U_fortsne);

gscatter(tsneresult(:,1),tsneresult(:,2),...
    extractAfter(extractBefore(U_names_fortsne,"_repall"),"iai"),...
    hexadecimalcolorcode2rgbtriplet(color),char(join(marker,'')))

end

