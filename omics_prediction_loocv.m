%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


%%%% Estimate proteomes from Raman spectra (LOOCV)


% For reproducibility
rng(1)

npermutations = 100000;

metric = 1;
% 1: PRESS
% 2: cosine distance
% 3: 1 - Pearson's correlation coefficient
% 4: L1 norm
% 5: median of relative errors


% LDA axes to use
meanU_axes = 1:4;

[npress, significance, proteinsmean_fgpercell_est, betas] = ...
    permutationtest(proteinsmean_fgpercell(:,:), meanU(:,meanU_axes), npermutations, metric);

disp("Estimation error: "+string(npress))
disp("p-value of permutation test: "+string(significance))


% condition to plot
pltcnd = 5;

figure;
hold on
plot(proteinsmean_fgpercell_est(pltcnd,:),proteinsmean_fgpercell(pltcnd,:),'.','Color',"#DAB24F")
plot([1e-8,1e+2],[1e-8,1e+2],'-','Color',"#DAB24F")
warning off MATLAB:Axes:NegativeDataInLogAxis % Estimated data contain negative values.
set(gca,'XScale','log')
set(gca,'YScale','log')
title(proteins_conditionnames(pltcnd))
xlabel("Estimated protein mass (fg/cell)")
ylabel("Measured protein mass (fg/cell)")

