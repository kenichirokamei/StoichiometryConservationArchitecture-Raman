%%% Copyright 2017-2023 Ken-ichiro F. Kamei %%%

% This code was created by revising a code in our previous paper (DOI: 10.1016/j.cels.2018.05.015).


%%%% Plot LDA Raman


% dimensions to plot
dims = [1,2];
if size(dims,2)~=2
    error('Specify two dimensions using a row vector.')
end

col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

figure;
plot_lda(U(:,2), dims, col, marr);
plot_lda_ellipse(U(:,2), dims, col);

h = findobj('Type','line');
set(h,'zdata',[])

title("LDA Raman")
xlabel(strcat('LDA',num2str(dims(1))))
ylabel(strcat('LDA',num2str(dims(2))))
legend(extractAfter(extractBefore(U(:,1),"_repall"),"iai"),'Interpreter','none')

