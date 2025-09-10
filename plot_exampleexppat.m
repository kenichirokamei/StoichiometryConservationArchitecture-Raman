%%% Copyright 2022-2023 Ken-ichiro F. Kamei %%%


function plot_exampleexppat(st_growthrate,st_ind,om,omse,om_description,pltprotein,xran,j,num)
    errorbar(st_growthrate,...
                om(st_ind,pltprotein(j)),omse(st_ind,pltprotein(j)),...
                '-o','Color',"#DAB24F",'MarkerSize',10)
    ylabel("Protein mass (fg/cell)")
    xlabel("Growth rate (1/min)")
    title(string(j+num)+": "+om_description(pltprotein(j),3))
    xlim(xran)
end

