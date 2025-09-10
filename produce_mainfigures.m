%%% Copyright 2023 Ken-ichiro F. Kamei %%%


% All codes were created by using MATLAB version 2019a Update 9 (9.6.0.1472908) 64-bit (maci64).


% All codes do not require any toolboxes 
% except for 'boxplot' in 'comp_cog' and 'visualize_ldaraman_tsne', 
% which use "Statistics and Machine Learning Toolbox" and are disabled by default.


% Description of the variables in the attached MAT-file:

    % Raman data --- This paper.
    % - 'rawraman'
    %     15-by-2 cell array; 1st column: condition name, 2nd column: Raman spectra.
    %     In each cell of the second column, each row represents a single-cell Raman spectrum. 

    % proteome data --- Extracted from supplementary tables 5 and 10 in A. Schmidt, et al., Nature Biotechnology 34, 104 (2016).
    % The following variables have the corresponding order of rows (15 conditions) and columns (2058 protein species).
    % - 'proteinsmean_fgpercell'
    %     15-by-2058 array; protein abundances in fg/cell.
    % - 'proteinscv'
    %     15-by-2058 array; protein abundance cv in %.
    % - 'proteins_conditionnames'
    %     15-by-1 array; condition names.
    % - 'proteins_description'
    %     3-by-2058 array; 1st row: Uniprot accession, 2nd row: explanation, 3rd row: gene name.
    % - 'COGclasses_index'
    %     2-by-4 cell array; 1st row: protein species indexes, 2nd row: COG class name.
    %     The protein species indexes are compatible with the order of the 2058 columns of the variables above. 

    % growthrate data --- Extracted from supplementary table 23 in A. Schmidt, et al., Nature Biotechnology 34, 104 (2016).
    % The following variables have the corresponding order of rows (15 conditions).
    % - 'growthrate'
    %     15-by-1 array; growth rates in 1/min.
    % - 'growthrate_conditionnames'
    %     15-by-1 array; condition names.

    % SCGs data --- Can also be created by executing 'extract_scgs'.
    % - 'SCGs_member_noden_original'
    %     5-by-3 cell array; 1st column: SCG name, 2nd column: member protein index, 3rd column: member protein explanation. 
    %     The protein indexes in the 2nd column are compatible with the proteome data mentioned above. 

    % The Raman data, the proteome data, and the growth data have the same order of the 15 conditions.


% Codes:

    % Applying LDA to Raman spectra
    % (Necessary to perform all analyses involving Raman data)
    raman_select_lda

        % Plotting LDA Raman (Fig. 2C–2F)
        plot_ldaraman
        %{
            % Fig. 2F (visualization of LDA Raman) uses functions 'tsne' and 'gscatter' in 
            % "Statistics and Machine Learning Toolbox."
            visualize_ldaraman_tsne
        %}

        % Predicting proteome data from Raman data, performing permutation test 
        % and observing coefficient matrix (Fig. 2H and 3A)
        omics_prediction_loocv
            scatterplot_coefficients % Dependent on 'omics_prediction_loocv'

        % Similarity of expression patterns between culture conditions for each COG class (Fig. 3B–3D)
        comp_cog % To create boxplots in Fig. 3B, 'boxplot' in "Statistics and Machine Learning Toolbox" is necessary (disabled by default.)
        constratio_examples
        
        % Extracting SCGs from proteome data with cosine similarity (Fig. 4B and 4C)
        extract_scgs
        
        % Structural similarity among the distribution of LDA Raman spectra, 
        % proteome structure determined by Raman-proteome transformation coefficients, 
        % and proteome structure determined by stoichiometry conservation (Fig. 5K, 5L, 6A–6D, S9A, S9F, S9G)
        link_ldaraman_csleomics

        % Stoichiometry conservation centrality and its proportionality to 
        % expression generality (Fig. 5A, 7A–7C)
        stoichioconscentrality_expressiongenerality
        
        % Dependence of low-dimensional correspondence between Raman spectra 
        % and proteomes on the number of conditions (Fig. S14 --- added during the revision at eLife)
        %{
            % This code is time-consuming and is therefore disabled by default. 
            subsample_cond
        %}

