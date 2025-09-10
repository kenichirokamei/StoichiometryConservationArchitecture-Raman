%%% Copyright 2021-2023 Ken-ichiro F. Kamei %%%


%%%% Structural similarity among the distribution of LDA Raman spectra, 
%%%% proteome structure determined by Raman-proteome transformation coefficients, 
%%%% and proteome structure determined by stoichiometry conservation (Fig. 5K, 5L, 6A-6D, S9A, S9F, S9G)


col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

selected_component_col = ["#A4538D",col(10),col(9),col(2),col(14)];
selected_component_marr = [".",marr(10),marr(9),marr(2),marr(14)];


% Load Raman
R = meanU;

% Load proteome
P = proteinsmean_fgpercell;
P_description = proteins_description';


ngroups = size(R,1);
nproteins = size(P,2);
labels = U(:,1);


R_E = [ones(ngroups,1) R];
S_R_E = diag([sqrt(ngroups);sqrt(diag(R'*R))]);
B_E = P'*(pinv(R_E))';

Phat = P./vecnorm(P,2,1);
deg = (Phat'*sum(Phat,2))';

[~,S_LE,V_sym] = svd(Phat./sqrt(deg),"econ");
    % For reproducibility, make the sign of the element with the largest absolute value in each
    % column positive.
    % (Signs of singular vectors that 'svd' function of Matlab outputs can be different 
    % across different machines and releases of Matlab.)
    [~,rowind_eachcol] = max(abs(V_sym),[],1);
    V_sym = V_sym .* sign(V_sym(sub2ind(size(V_sym),rowind_eachcol,1:size(V_sym,2))));
V_rw = (1./sqrt(deg')).*V_sym;
B_E_est = (vecnorm(P,2,1).*sqrt(deg))'.*V_sym;

theta = pinv(B_E_est.*diag(S_LE)')*B_E.*diag(S_R_E)';


figure;
hold on
scatter(B_E(:,1)*S_R_E(1,1),B_E_est(:,1)*S_LE(1,1),'o','MarkerEdgeColor',"#DAB24F")
plot(xlim,xlim,'-.','Color',"#AAAAAA")
xlabel("$\sqrt{m}\ B_{E}[:,1]$",'Interpreter','latex')
ylabel("$B_{E}^{\mathrm{est}}[:,1]$",'Interpreter','latex')

figure;
b = bar3(abs(theta));
for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    b(i).FaceAlpha = 0.5;
end
colorbar
caxis([-1 1])
colormap cool
title("$|\Theta|$",'Interpreter','latex')
xlabel("column")
ylabel("row")
zlim([-1 1])
view(25,30)


xax = 1;
yax = 2;
zax = 4;
figure;
scatter3(V_rw(:,xax+1),V_rw(:,yax+1),V_rw(:,zax+1),[],deg,'.','SizeData',30)
xlabel("csLE"+string(xax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
ylabel("csLE"+string(yax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
zlabel("csLE"+string(zax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(zax+1)+"]$)",'Interpreter','latex')
title(["Proteome","stoichiometric balance"])
c = colorbar;
c.Label.String = "Stoichiometry conservation centrality";
colormap(gca,"parula")
grid on
view(30,10)
set(gca,'YDir','reverse')

xax = 1;
yax = 2;
zax = 4;
figure('Position', [0 0 1000 420]);
s1 = subplot(1,2,1);
    hold on
    scatter3(V_rw(:,xax+1),V_rw(:,yax+1),V_rw(:,zax+1),...
        'filled','MarkerFaceColor',"#777777",'MarkerFaceAlpha',0.2,'SizeData',5)
    p = [];
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter3(V_rw(ind,xax+1),V_rw(ind,yax+1),V_rw(ind,zax+1),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1});
    end
    legend(p,'Location','northwest');
    xlabel("csLE"+string(xax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
    ylabel("csLE"+string(yax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
    zlabel("csLE"+string(zax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(zax+1)+"]$)",'Interpreter','latex')
    grid on
    view(30,10)
    set(gca,'YDir','reverse')
s2 = subplot(1,2,2);
    copyobj(get(s1,'children'),s2)
    labeltmp = get(s1,'xlabel');
    xlabel(labeltmp.String,'Interpreter','latex')
    labeltmp = get(s1,'ylabel');
    ylabel(labeltmp.String,'Interpreter','latex')
    labeltmp = get(s1,'zlabel');
    zlabel(labeltmp.String,'Interpreter','latex')
    grid on
    view(80,20)
    set(gca,'YDir','reverse')
sgtitle(["Proteome","stoichiometric balance"])

xax = 1;
yax = 2;
zax = 3;
figure('Position', [0 0 1000 420]);
s1 = subplot(1,2,1);
    hold on
    scatter3(B_E(:,xax+1)./B_E(:,1),B_E(:,yax+1)./B_E(:,1),B_E(:,zax+1)./B_E(:,1),...
        'filled','MarkerFaceColor',"#777777",'MarkerFaceAlpha',0.2,'SizeData',5)
    p = [];
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter3(B_E(ind,xax+1)./B_E(ind,1),B_E(ind,yax+1)./B_E(ind,1),B_E(ind,zax+1)./B_E(ind,1),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1});
    end
    legend(p,'Location','northwest');
    xlabel("$\beta^{\mathrm{LDA}"+xax+"}$ ($B_E^{\mathrm{norm}}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
    ylabel("$\beta^{\mathrm{LDA}"+yax+"}$ ($B_E^{\mathrm{norm}}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
    zlabel("$\beta^{\mathrm{LDA}"+zax+"}$ ($B_E^{\mathrm{norm}}[:,"+string(zax+1)+"]$)",'Interpreter','latex')
    grid on
    view(20,25)
s2 = subplot(1,2,2);
    copyobj(get(s1,'children'),s2)
    labeltmp = get(s1,'xlabel');
    xlabel(labeltmp.String,'Interpreter','latex')
    labeltmp = get(s1,'ylabel');
    ylabel(labeltmp.String,'Interpreter','latex')
    labeltmp = get(s1,'zlabel');
    zlabel(labeltmp.String,'Interpreter','latex')
    grid on
    view(60,20)
sgtitle(["Proteome","Raman-proteome coef."])

xax = 1;
yax = 2;
figure('Position', [0 0 1680 1000]);
subplot(2,3,1);
    hold on
    p = [];
    for i=1:ngroups
        p(i) = plot(R_E(i,xax+1),R_E(i,yax+1),...
            marr(i),'Color',col(i),'DisplayName',string(extractBetween(labels{i},'iai','_')));
    end
    xline(0,'--','Color','#AAAAAA');
    yline(0,'--','Color','#AAAAAA');
    legend(p,'Location','south');
    xlabel("LDA"+string(xax))
    ylabel("LDA"+string(yax))
    title(["Raman","LDA"])
subplot(2,3,2);
    hold on
    scatter(B_E(:,xax+1)./B_E(:,1),B_E(:,yax+1)./B_E(:,1),...
        'filled','MarkerFaceColor',"#AAAAAA",'MarkerFaceAlpha',0.5,'SizeData',7)
    p = [];
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter(B_E(ind,xax+1)./B_E(ind,1),B_E(ind,yax+1)./B_E(ind,1),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1});
    end
    xline(0,'--','Color','#AAAAAA');
    yline(0,'--','Color','#AAAAAA');
    legend(p,'Location','south');
    xlabel("$\beta^{\mathrm{LDA}"+xax+"}$ ($B_E^{\mathrm{norm}}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
    ylabel("$\beta^{\mathrm{LDA}"+yax+"}$ ($B_E^{\mathrm{norm}}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
    title(["Proteome","Raman-proteome coef."])
subplot(2,3,3);
    hold on
    scatter(V_rw(:,xax+1),V_rw(:,yax+1),...
        'filled','MarkerFaceColor',"#AAAAAA",'MarkerFaceAlpha',0.5,'SizeData',7)
    p = [];
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter(V_rw(ind,xax+1),V_rw(ind,yax+1),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1});
    end
    xline(0,'--','Color','#AAAAAA');
    yline(0,'--','Color','#AAAAAA');
    legend(p,'Location','south');
    xlabel("csLE"+string(xax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
    ylabel("csLE"+string(yax)+" ($\tilde{V}_\mathrm{rw}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
    title(["Proteome","stoichiometric balance"])
    set(gca,'YDir','reverse')
subplot(2,3,6);
    hold on
    scatter(B_E_est(:,xax+1)./B_E_est(:,1),B_E_est(:,yax+1)./B_E_est(:,1),...
        'filled','MarkerFaceColor',"#AAAAAA",'MarkerFaceAlpha',0.5,'SizeData',7)
    p = [];
    for i=1:size(SCGs_member_noden_original,1)
        ind = SCGs_member_noden_original{i,2};
        p(end+1) = scatter(B_E_est(ind,xax+1)./B_E_est(ind,1),B_E_est(ind,yax+1)./B_E_est(ind,1),...
            selected_component_marr(mod(i-1,length(selected_component_marr))+1),...
            'MarkerEdgeColor',selected_component_col(mod(i-1,length(selected_component_col))+1),...
            'DisplayName',SCGs_member_noden_original{i,1});
    end
    xline(0,'--','Color','#AAAAAA');
    yline(0,'--','Color','#AAAAAA');
    legend(p,'Location','south');
    xlabel("csLE"+string(xax)+" ($B_E^{\mathrm{est,norm}}[:,"+string(xax+1)+"]$)",'Interpreter','latex')
    ylabel("csLE"+string(yax)+" ($B_E^{\mathrm{est,norm}}[:,"+string(yax+1)+"]$)",'Interpreter','latex')
    title(["On the basis of the mathematical analysis,",...
        "this csLE plot uses $B_E^{\mathrm{est,norm}} = \left(\sum_i d_i\right)^{1/2}\tilde{V}_\mathrm{rw} = \left(\sum_i d_i\right)^{1/2}D^{-1/2}\tilde{V}_\mathrm{sym}$, where each column of $\tilde{V}_\mathrm{sym}$ is normalized.",...
        "See supplementary materials for details."],'interpreter','latex')
    set(gca,'YDir','reverse')

