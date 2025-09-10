%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


function plot_lda_ellipse(U, dims, color)

m = {'-'};
set(gca,'LineStyleOrder',m,'ColorOrder',hexadecimalcolorcode2rgbtriplet(color),'ColorOrderIndex',1);
hold all

ngroups = size(U,1);
for i=1:ngroups
    cova = cov(U{i}(:,dims));
    me = mean(U{i}(:,dims));

    error_ellipse(cova,me,'conf',0.95);
    % Downloaded from https://jp.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse
end

end

