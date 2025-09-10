%%% Copyright 2021-2023 Ken-ichiro F. Kamei %%%


function plot_scgs_overview(G_abovethresh,COGclasses_index_wn,coloredCOGclass)
    
    if length(coloredCOGclass)~=1
        error('Specify a single COG class.')
    end

    COGclasses = table2array(cell2table(COGclasses_index_wn(2,:)));
    COGclasses_index = COGclasses_index_wn(1,:);

    coloredCOGclasslabel = [lower(COGclasses(coloredCOGclass))+" COG class","Other"];
    
    colormapvec = hexadecimalcolorcode2rgbtriplet(["#777777","#AAAAAA"]);
    colorbar_pos = linspace(1,3,size(colormapvec,1)*2+1);
    G_abovethresh.Nodes.NodeColorIndexes = repmat(colorbar_pos(end),size(G_abovethresh.Nodes,1),1);
    
    for i=1:length(COGclasses)
        for j=1:size(COGclasses_index{i},2)
            G_abovethresh.Nodes.COGclasses(G_abovethresh.Nodes.Noden_original==COGclasses_index{i}(j)) = COGclasses(i);
            if i==coloredCOGclass
                G_abovethresh.Nodes.NodeColorIndexes(G_abovethresh.Nodes.Noden_original==COGclasses_index{i}(j)) = colorbar_pos(1);
            end
        end
    end
    
    figure;
    h2 = plot(G_abovethresh,...
        'NodeCData',G_abovethresh.Nodes.NodeColorIndexes,...
        'MarkerSize',3);
    colormap(colormapvec)
    colorbar('eastoutside',...
        'Ticks',colorbar_pos(2:2:end),...
        'TickLabels',coloredCOGclasslabel,...
        'TickLength',0,...
        'Location','southoutside')
    h2.EdgeColor = "#555555";
    layout(h2,'force','UseGravity',true)

end

