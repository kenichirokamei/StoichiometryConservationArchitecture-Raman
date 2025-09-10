%%% Copyright 2022-2023 Ken-ichiro F. Kamei %%%


%%%% Examples of stoichiometry-conserving proteins in the ISP COG class and examples of non-ISP COG class proteins


col = ["#ACC255","#DAB24F",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C","#A4538D","#ACC255","#3E90BA",...
    "#DAB24F","#C35866","#7C80AA","#37A372","#E5D64C"];
marr = ["*","x","s","d","^","v",">","<","p","h","+",...
    "*","x","s","d"];

nreplicates = 3; % According to A. Schmidt, et al., Nature Biotechnology 34, 104 (2016).
proteinsse = proteinscv.*proteinsmean_fgpercell/100/sqrt(nreplicates);


%%% some ISP COG class proteins in the homeostatic core

% 375: rplF (50S ribosomal subunit protein L6)
% 866: rimM (ribosome maturation factor RimM)
% 42: tsf (protein chain elongation factor EF-Ts)
% 1304: rluB (23S rRNA pseudouridine2605 synthase)
% 30: aspS (aspartate—tRNA ligase)
% 94: rho (transcription termination factor Rho)
% 73: polA (DNA polymerase I)
% (Description are extracted from EcoCyc on August 15th, 2022.)

selprot = [375,866,42,1304,30,94,73];
selprot_desc = ["50S ribosomal subunit protein L6",...
    "ribosome maturation factor RimM",...
    "protein chain elongation factor EF-Ts",...
    "23S rRNA pseudouridine2605 synthase",...
    "aspartate-tRNA ligase",...
    "transcription termination factor Rho",...
    "DNA polymerase I"]; 

figure;
hold on
stdprot = 1; % rplF
xmin = 0;
xmax = 1.7;
p = [];
for i=setdiff(1:length(selprot),stdprot)
    p = [p,...
        errorbar(proteinsmean_fgpercell(:,selprot(stdprot)),... % x
        proteinsmean_fgpercell(:,selprot(i)),... % y
        proteinsse(:,selprot(i)),... % yneg
        proteinsse(:,selprot(i)),... % ypos
        proteinsse(:,selprot(stdprot)),... % xneg
        proteinsse(:,selprot(stdprot)),... % xpos
        marr(mod(i-1-1,length(marr))+1),'Color',col(mod(i+8-1,length(col))+1),'MarkerSize',10,...
        'DisplayName',proteins_description(3,selprot(i))+" ("+selprot_desc(i)+")")];

    reg = NaN(2,1);
    % least-squares regression (intersept: 0)
    reg(2,:) =...
        sum(proteinsmean_fgpercell(:,selprot(stdprot)).*proteinsmean_fgpercell(:,selprot(i)))/...
        sum(proteinsmean_fgpercell(:,selprot(stdprot)).^2);
    reg(1,:) = 0;
    plot([xmin,xmax],[reg(1,:)+xmin*reg(2,:),reg(1,:)+xmax*reg(2,:)],'-','Color',col(mod(i+8-1,length(col))+1))
end
xlim([0 xmax])
ylim([0 2.5])
lgd = legend(p,'Location','southoutside');
title(lgd,'Compared gene')
xlabel(["Abundance (fg/cell)",proteins_description(3,selprot(stdprot))+" ("+selprot_desc(stdprot)+")"])
ylabel(["Abundance (fg/cell)","compared gene"])


%%% some non-ISP COG class proteins whose abundance ratios are not constant

% 375: rplF (50S ribosomal subunit protein L6)
% 350: crp (DNA-binding transcriptional dual regulator CRP)
% 249: fbaB (fructose-bisphosphate aldolase class I, in glycolysis)
% 72: acnA (aconitate hydratase A, in TCA cycle)
% 388: aroA (EPSP synthase, in shikimate pathway)
% 306: dapB (DHDPR, in diaminopimelate pathway)
% 4: purL (phosphoribosylformylglycinamide synthetase, in purine de novo biosynthesis pathway)
% (Description are extracted from EcoCyc on August 22nd, 2022.)

selprot = [375,350,249,72,388,306,4];
selprot_desc = ["50S ribosomal subunit protein L6",...
    "DNA-binding transcriptional dual regulator CRP",...
    "in glycolysis",...
    "in TCA cycle",...
    "in shikimate pathway",...
    "in diaminopimelate pathway",...
    "in purine de novo biosynthesis pathway"];

figure;
hold on
stdprot = 1; % rplF
xmin = 0;
xmax = 1.7;
p = [];
[~,sortind] = sort(proteinsmean_fgpercell(:,selprot(stdprot)));
for i=setdiff(1:length(selprot),stdprot)
    p = [p, errorbar(proteinsmean_fgpercell(sortind,selprot(stdprot)),... % x
        proteinsmean_fgpercell(sortind,selprot(i)),... % y
        proteinsse(sortind,selprot(i)),... % yneg
        proteinsse(sortind,selprot(i)),... % ypos
        proteinsse(sortind,selprot(stdprot)),... % xneg
        proteinsse(sortind,selprot(stdprot)),... % xpos
        marr(mod(i-1-1,length(marr))+1)+"-.",'Color',col(mod(i+8-1,length(col))+1),'MarkerSize',10,...
        'DisplayName',proteins_description(3,selprot(i))+" ("+selprot_desc(i)+")")];
end
xlim([0 xmax])
ylim([0 0.65])
lgd = legend(p,'Location','southoutside');
title(lgd,'Compared gene')
xlabel(["Abundance (fg/cell)",proteins_description(3,selprot(stdprot))+" ("+selprot_desc(stdprot)+")"])
ylabel(["Abundance (fg/cell)","compared gene"])

