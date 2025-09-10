%%% Copyright 2020-2023 Ken-ichiro F. Kamei %%%


%%%% analysis of betas


whichCOGclass = 2; % Information storage and processing (ISP) COG class

index2 = COGclasses_index{1,whichCOGclass};
selected_COGclass_name = COGclasses_index{2,whichCOGclass};

nconds = size(proteinsmean_fgpercell,1);
nproteins = size(proteinsmean_fgpercell,2);

% make "B" as the average of betas
B = zeros(size(proteinsmean_fgpercell,2),length(meanU_axes)+1);
for prot=1:nproteins
    for i=1:length(meanU_axes)+1
        for cond=1:nconds
            B(prot,i) = B(prot,i) + betas{cond,prot}(i);
        end
    end
end
B = B/nconds;


figure;
xx = 0;
tmp = 0;
for yy=meanU_axes
    tmp = tmp + 1;
    
    subplot(2,2,tmp)
    hold on
        p = [];
        p = [p, plot(B(:,xx+1),B(:,yy+1),'.','Color',"#3E90BA",'DisplayName','All')];
        p = [p, plot(B(index2,xx+1),B(index2,yy+1),'o','Color',"#DAB24F",...
            'DisplayName',lower(selected_COGclass_name)+" COG class")];

            % least-squares regression (intercept: 0)
            sl = B(index2,xx+1)\B(index2,yy+1);
            reg_xlim = [0 20];
            p = [p, plot(reg_xlim,sl*reg_xlim,'-','Color',"#DAB24F")];

        plot([-10 20],[0 0],'--','Color','k')
        plot([0,0],[-10 20],'--','Color','k')
        
        xlim([floor(min(B(:,xx+1))) reg_xlim(2)])
        ylim([-10 20])

        xlabel("Constants")
        ylabel("Coefficients for LDA"+string(yy))

        if tmp==1
            legend(p(1:2))
        end
end

