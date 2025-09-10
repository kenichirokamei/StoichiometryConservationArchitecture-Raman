%%% Copyright 2018-2023 Ken-ichiro F. Kamei %%%


%%%% Visalize LDA Raman by t-SNE


rng('shuffle')

col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

figure;
plot_lda_tsne(U, col, marr)

title("t-SNE of LDA Raman")

