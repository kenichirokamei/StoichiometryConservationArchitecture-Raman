%%% Copyright 2022-2023 Ken-ichiro F. Kamei %%%


%%%% Stoichiometry conservation centrality and its proportionality to expression generality


col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

selected_component_col = ["#A4538D",col(10),col(9),col(2),col(14)];
selected_component_marr = [".",marr(10),marr(9),marr(2),marr(14)];



P = proteinsmean_fgpercell;
P_description = proteins_description';

nreplicates = 3; % According to A. Schmidt, et al., Nature Biotechnology 34, 104 (2016).
Pse = proteinscv.*P/100/sqrt(nreplicates);


%%% L1 norm/L2 norm vs. degree

% expression generality
l1overl2 = vecnorm(P,1,1)./vecnorm(P,2,1);

% stoichiometry conservation centrality 
% (In our case, it is equal to the eigenvector corresponding to the largest 
% eigenvalue of a "normalized" adjacency matrix. See supplementary materials.)
Phat = P./vecnorm(P,2,1);
deg = (Phat'*sum(Phat,2))';

% example proteins
l1overl2_low_pltprotein = [2039,1767,1769]; % low L1/L2 from top to bottom
l1overl2_mid_pltprotein = [1692,1489,1377]; % middle L1/L2 from top to bottom
l1overl2_mid2_pltprotein = [554,1982,79]; % higher middle L1/L2 from top to bottom
l1overl2_high_pltprotein = [505,1216,1199]; % high L1/L2 from top to bottom

% growth rates
[sorted_growthrate,sorted_ind] = sort(growthrate);

figure('Position',[0 0 1700 1000]);
subplot(3,7,[1,2,3,8,9,10,15,16,17])
    hold on
    scatter(l1overl2,deg,'filled','MarkerFaceColor',"#AAAAAA",'MarkerFaceAlpha',0.5,'SizeData',10)
    p = [];
    plot([0 ceil(sqrt(size(P,1)))],[0 ceil(sqrt(size(P,1)))]*sqrt(sum(deg))/sqrt(size(P,1)),'-','Color',"#AAAAAA",'LineWidth',1)
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter(l1overl2(ind),deg(ind),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1},'SizeData',80);
    end
    txtposx = -0.1;
    txtposy = 50;
    scatter(l1overl2(l1overl2_low_pltprotein),deg(l1overl2_low_pltprotein),'o','SizeData',100,'MarkerEdgeColor','r')
    text(l1overl2(l1overl2_low_pltprotein)+txtposx,deg(l1overl2_low_pltprotein)+txtposy,string(1:3)'+": "+P_description(l1overl2_low_pltprotein,3),'FontSize',15)
    scatter(l1overl2(l1overl2_mid_pltprotein),deg(l1overl2_mid_pltprotein),'o','SizeData',100,'MarkerEdgeColor','r')
    text(l1overl2(l1overl2_mid_pltprotein)+txtposx,deg(l1overl2_mid_pltprotein)+txtposy,string(4:6)'+": "+P_description(l1overl2_mid_pltprotein,3),'FontSize',15)
    scatter(l1overl2(l1overl2_mid2_pltprotein),deg(l1overl2_mid2_pltprotein),'o','SizeData',100,'MarkerEdgeColor','r')
    text(l1overl2(l1overl2_mid2_pltprotein)+txtposx,deg(l1overl2_mid2_pltprotein)+txtposy,string(7:9)'+": "+P_description(l1overl2_mid2_pltprotein,3),'FontSize',15)
    scatter(l1overl2(l1overl2_high_pltprotein),deg(l1overl2_high_pltprotein),'o','SizeData',100,'MarkerEdgeColor','r')
    text(l1overl2(l1overl2_high_pltprotein)+txtposx,deg(l1overl2_high_pltprotein)+txtposy,string(10:12)'+": "+P_description(l1overl2_high_pltprotein,3),'FontSize',15)
    xline(1,'--','Color',"#AAAAAA");
    xline(sqrt(size(P,1)),'--','Color',"#AAAAAA");
    yline(size(P,2),'--','Color',"#AAAAAA");
    xrange = xlim;
    yrange = ylim;
    lgd = legend(p,'Location','southeast');
    lgd.NumColumns = 1;
    title(lgd,"Expressed Mainly in ...")
    xlim([0 xrange(2)])
    ylim([0 2100])
    xlabel(["Expression generality", "(L1 norm / L2 norm of expression vector)"])
    ylabel(["Stoichiometry conservation centrality","(sum of cosine similarity)"])

xrange = [-0.005 0.035];
for i=1:length(l1overl2_low_pltprotein)
    subplot(3,7,(i-1)*7+4)
        plot_exampleexppat(sorted_growthrate,sorted_ind,P,Pse,P_description,l1overl2_low_pltprotein,xrange,i,0)
end
for i=1:length(l1overl2_mid_pltprotein)
    subplot(3,7,(i-1)*7+5)
        plot_exampleexppat(sorted_growthrate,sorted_ind,P,Pse,P_description,l1overl2_mid_pltprotein,xrange,i,3)
end
for i=1:length(l1overl2_mid2_pltprotein)
    subplot(3,7,(i-1)*7+6)
        plot_exampleexppat(sorted_growthrate,sorted_ind,P,Pse,P_description,l1overl2_mid2_pltprotein,xrange,i,6)
end
for i=1:length(l1overl2_high_pltprotein)
    subplot(3,7,(i-1)*7+7)
        plot_exampleexppat(sorted_growthrate,sorted_ind,P,Pse,P_description,l1overl2_high_pltprotein,xrange,i,9)
end


%%% stoichiometry conservation centrality distribution

binwidth = 100;

figure;
hold on
    % all
    histogram(deg,'BinWidth',binwidth,'Normalization','probability',...
        'FaceColor',"#AAAAAA",'DisplayName',"All")
    % homeostatic core
    i=1;
    ind_hc = SCGs_member_noden_original{i,2};
    histogram(deg(ind_hc),'BinWidth',binwidth,'Normalization','probability',...
        'FaceColor',selected_component_col(i),'DisplayName',SCGs_member_noden_original{i,1})
    % other SCGs
    i=2:size(SCGs_member_noden_original,1);
    ind_otherSCGs = cell2mat(SCGs_member_noden_original(i,2));
    histogram(deg(ind_otherSCGs),'BinWidth',binwidth,'Normalization','probability',...
        'FaceColor',selected_component_col(2),...
        'DisplayName',join(table2array(cell2table(SCGs_member_noden_original(i,1))),", "))
legend('show','Location','northwest')
xlabel("Stoichiometry conservation centrality")
ylabel("Relative frequency")
ylim([0 1])

