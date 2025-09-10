%%% Copyright 2025 Ken-ichiro F. Kamei %%%


function [pcc] = pcc_squaremat(P)
    
    Psize = size(P,1);
    Prowcolnum = 1:Psize;
    
    samplesize = sum(P,'all');
    xsum = Prowcolnum*sum(P,2);
    xsqsum = Prowcolnum.^2*sum(P,2);
    ysum = sum(P,1)*Prowcolnum';
    ysqsum = sum(P,1)*(Prowcolnum').^2;
    xysum = Prowcolnum*P*Prowcolnum';
    
    pcc = (samplesize*xysum - xsum*ysum) / (sqrt(samplesize*xsqsum - xsum^2) * sqrt(samplesize*ysqsum - ysum^2));

end

