%%% Copyright 2025 Ken-ichiro F. Kamei %%%


%%%% Dependence of low-dimensional correspondence between Raman spectra and proteomes 
%%%% on the number of conditions --- added during the revision at eLife


col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];


% For reproducibility
rng(1)

randorthog_nsamplings = 10000;

significances = cell(ngroups,1);
for nconds=2:ngroups
    subsamples = nchoosek(1:ngroups,nconds);
    nsubsamples = nchoosek(ngroups,nconds);
    
    pcc_randorthog = NaN(randorthog_nsamplings,1);
    for j=1:randorthog_nsamplings
        [QQ,RR] = qr(randn(nconds,nconds));
        randorthog = QQ * diag(sign(diag(RR)));
        pcc_randorthog(j) = pcc_squaremat(randorthog.^2);
    end
    
    pcc_theta = NaN(nsubsamples,1);
    significance_tmp = NaN(nsubsamples,1);
    for i=1:nsubsamples
        
         %%% Conduct LDA
        U = lda(selraman(subsamples(i,:),2), selraman(subsamples(i,:),1), 98);

        % calculate mean of each group
        meanU = NaN(nconds, nconds-1);
        for j=1:nconds
            meanU(j,:) = mean(U{j,2}, 1);
        end
        
        % Load Raman
        R = meanU;
        % Load proteome
        P = proteinsmean_fgpercell(subsamples(i,:),:);

        R_E = [ones(nconds,1) R];
        S_R_E = diag([sqrt(nconds);sqrt(diag(R'*R))]);
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
        B_E_est = (vecnorm(P,2,1).*sqrt(deg))'.*V_sym;

        theta = pinv(B_E_est.*diag(S_LE)')*B_E.*diag(S_R_E)';

        pcc_theta(i) = pcc_squaremat(theta.^2);
        significance_tmp(i) = sum(pcc_randorthog>pcc_theta(i))/randorthog_nsamplings;
    end
    significances{nconds} = significance_tmp;
end


figure;
hold on
p = [];
for nconds=2:ngroups
    tmp = significances{nconds};
    tmp(significances{nconds}==0) = 1/randorthog_nsamplings/10;
    plot(ones(nchoosek(ngroups,nconds),1)*nconds,tmp,'s','Color',col(1))
    tmp = median(significances{nconds});
    if tmp==0
        tmp = 1/randorthog_nsamplings/10;
    end
    if nconds==2
        p(end+1) = plot([nconds-0.3 nconds+0.3], repmat(tmp,1,2),'-','LineWidth',1,'Color',[0.1 0.1 0.1],...
                    'DisplayName',"Median");
    else
        plot([nconds-0.3 nconds+0.3], repmat(tmp,1,2),'-','LineWidth',1,'Color',[0.1 0.1 0.1])
    end
end
detectionlim = 1/randorthog_nsamplings;
p(end+1) = yline(detectionlim,'--','LineWidth',1,'Color',col(10),'DisplayName',"Detection limit");
xlabel("#conditions used")
ylabel(["Probability of accidentally obtaining higher level of","low-dimensional correspondence than", ...
    "that of experimental data"])
xticks(2:ngroups)
xticklabels(string(2:ngroups))
yticks(10.^(-log10(randorthog_nsamplings)-1:-1))
yticklabels(["No occurrence" string([10.^(-log10(randorthog_nsamplings):-1)])])
legend(p,'Location','northeast')
xlim([2-0.5 ngroups+0.5])
ylim([1/randorthog_nsamplings/10/2 max(cell2mat(significances))*2])
set(gca,'YScale','log')
set(gca,'YMinorTick','off')
set(gca,'FontSize',15)

