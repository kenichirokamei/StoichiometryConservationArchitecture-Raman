%%% Copyright 2020-2023 Ken-ichiro F. Kamei %%%


%%%% Similarity of expression patterns between culture conditions for each COG class


% Boxplots in Fig. 3B uses a function 'boxplot' in "Statistics and Machine Learning Toolbox."
% Boxplots are created when 'boxp' is "on".
boxp = "off";


nconds = size(proteinsmean_fgpercell,1);

COGclasses_proteinsmean_fgpercell = cell(1,4);

for j=1:size(COGclasses_index,2)
    index2 = COGclasses_index{1,j};
    COGclasses_proteinsmean_fgpercell{1,j} = proteinsmean_fgpercell(:,index2);
end

COGClasses_names = table2array(cell2table(COGclasses_index(2,1:3)));

corr_coef = NaN(nchoosek(nconds,2),length(COGClasses_names));
corr_coef_log = NaN(nchoosek(nconds,2),length(COGClasses_names));

for i=1:length(COGClasses_names)
    corr_coef(:,i) = corr_pearson(COGclasses_proteinsmean_fgpercell{i}','vector');
    corr_coef_log(:,i) = corr_pearson(log10(COGclasses_proteinsmean_fgpercell{i}'),'vector');
end

figure;
hold on
if boxp=="on"
    boxplot(corr_coef(:,1:3),'Widths',0.2,'Colors',[0.3,0.3,0.3],'Notch','on')
end
for i=1:length(COGClasses_names)
    plot(i+0.2,corr_coef(:,i),'x','Color',col(i+7))
end
xlabel("COG class")
ylabel(["Pearson's correlation coefficients","between conditions"])
xlim([0 length(COGClasses_names)+1])
xticks(1:length(COGClasses_names))
xticklabels(lower(COGClasses_names))
xtickangle(45)

figure;
hold on
if boxp=="on"
    boxplot(corr_coef_log(:,1:3),'Widths',0.2,'Colors',[0.3,0.3,0.3],'Notch','on')
end
for i=1:length(COGClasses_names)
    plot(i+0.2,corr_coef_log(:,i),'x','Color',col(i+7))
end
xlabel("COG class")
ylabel(["Pearson's correlation coefficients","between conditions (log10 abundance)"])
xlim([0 length(COGClasses_names)+1])
xticks(1:length(COGClasses_names))
xticklabels(lower(COGClasses_names))
xtickangle(45)

