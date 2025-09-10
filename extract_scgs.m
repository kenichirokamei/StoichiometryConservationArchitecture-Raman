%%% Copyright 2020-2023 Ken-ichiro F. Kamei %%%


% Extract SCGs from proteome data with cosine similarity


col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

nreplicates = 3; % According to A. Schmidt, et al., Nature Biotechnology 34, 104 (2016).
proteinsse = proteinscv.*proteinsmean_fgpercell/100/sqrt(nreplicates);

% adjacency matrix (The diagonal elements are ones for simplicity.)
A = (proteinsmean_fgpercell./vecnorm(proteinsmean_fgpercell,2,1))'*...
    proteinsmean_fgpercell./vecnorm(proteinsmean_fgpercell,2,1);


% make a proteome graph based on stoichiometry conservation relations
G = graph(A,'upper');
nproteins = size(proteinsmean_fgpercell,2);
G = rmedge(G,1:nproteins,1:nproteins); % remove self-loops
G.Nodes.Noden_original = [1:nproteins]'; % record node numbers.

thresh = 0.995; % cosine similarity threshold for extracting SCGs
G_abovethresh = rmedge(G,find(G.Edges.Weight<thresh)); % remove edges of which weight is below the threshold and make a new graph
G_abovethresh = rmnode(G_abovethresh,find(degree(G_abovethresh)==0)); % remove isolated nodes 

component_index = conncomp(G_abovethresh)'; % detect components
    % For reproducibility, redefine component numbers.
    % (The documentation of 'conncomp' function of Matlab does not guarantee
    % that the order of the component numbering that the function outputs is 
    % the same across different machines and releases of Matlab.
    % The following code reproduces the output by Matlab 2019a Update 9 (9.6.0.1472908) 64-bit (maci64)
    % on Macbook Pro 13-inch, M2, 2022.)
    component_index_new = component_index;
    ncomp = length(unique(component_index));
    comp_firstpos = NaN(1,ncomp);
    for i=1:ncomp
        comp_firstpos(i) = find(component_index==i,1);
    end
    sorted_comp_firstpos = sort(comp_firstpos);
    for i=1:ncomp
        component_index_new(component_index==component_index(sorted_comp_firstpos(i))) = i;
    end
G_abovethresh.Nodes.Componentn = component_index_new;

whichCOGclass = 2; % Information storage and processing (ISP) COG class
plot_scgs_overview(G_abovethresh,COGclasses_index,whichCOGclass);

selected_component = [1,12,42,22,38];
% 'selected_component' is made for 'component_index_new' ('thresh': 0.995).
    % component 1: homeostatic core (SCG1)
    % component 12: mainly expressed under LB (SCG2)
    % component 42: mainly expressed under GlycerolAA (SCG3)
    % compoennt 22: mainly expressed under Fructose (SCG4)
    % component 38: mainly expressed under Stationray3days (SCG5)

SCGs_member_noden_original = cell(length(selected_component),3);
for i=1:length(selected_component)
    SCGs_member_noden_original{i,1} = "SCG"+string(i);
    SCGs_member_noden_original{i,2} =...
        G_abovethresh.Nodes.Noden_original(G_abovethresh.Nodes.Componentn==selected_component(i));
    SCGs_member_noden_original{i,3} =...
        proteins_description([3,2],SCGs_member_noden_original{i,2})';
    SCGs_member_noden_original{i,3}(:,2) =...
        extractBefore(SCGs_member_noden_original{i,3}(:,2)," OS=");
end

selected_component_col = ["#A4538D",col(10),col(9),col(2),col(14)];
selected_component_mar = [".",marr(10),marr(9),marr(2),marr(14)]+"-";

[sorted_growthrate,sorted_ind] = sort(growthrate);
figure;
for i=1:length(selected_component)
    subplot(2,ceil(length(selected_component)/2),i)
    hold on
    for j=1:length(SCGs_member_noden_original{i,2})
        errorbar(sorted_growthrate,...
            proteinsmean_fgpercell(sorted_ind,SCGs_member_noden_original{i,2}(j)),...
            proteinsse(sorted_ind,SCGs_member_noden_original{i,2}(j)),...
            selected_component_mar(mod(i-1,length(selected_component_mar))+1),...
            'Color',selected_component_col(mod(i-1,length(selected_component_col))+1))
    end
    title(SCGs_member_noden_original{i,1})
    xlim([-0.001 0.035])
    xlabel("Growth rate (1/min)")
    ylabel("Protein mass (fg/cell)")
end

