function [TP,FN,FP,TN,n_GT,n_BW,A,B,C,D] = ComputeContingencyTable(GT,BW)
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    GT = logical(GT);
    BW = logical(BW);
    TP = 0;
    FN = 0;
    FP = 0;
    TN = 0;
    n_GT = 0;
    n_BW = 0;
    [x,y] = size(GT);
    for i=1:x
        for j=1:y
            if(GT(i,j) == 1)
                n_GT = n_GT + 1;
            end
            if(BW(i,j) == 1)
                n_BW = n_BW + 1;
            end
            if(GT(i,j) == 1 && BW(i,j) == 1)
                TP = TP + 1;
            end
            if(GT(i,j) == 1 && BW(i,j) == 0)
                FN = FN + 1;
            end
            if(GT(i,j) == 0 && BW(i,j) == 1)
                FP = FP + 1;
            end
            if(GT(i,j) == 0 && BW(i,j) == 0)
                TN = TN + 1;
            end
        end
    end
    coltot1 = TN + FP;
    coltot2 = FN + TP;
    rowtot1 = TN + FN;
    rowtot2 = TP + FP;
    nis = rowtot1^2 + rowtot2^2;
    njs = coltot1^2 + coltot2^2;
    s = TP^2 + TN^2 + FP^2 + FN^2;
    A = (Binomial(TN) + Binomial(FN) + Binomial(FP) + Binomial(TP) ) / 2;
    B = ( njs - s ) / 2;
    C = ( nis - s ) / 2;
    D = ( (x*y)^2 + s - nis - njs ) / 2;
end
