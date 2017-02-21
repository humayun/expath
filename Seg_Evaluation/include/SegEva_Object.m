function [S_TP,S_FN,S_FP,S_TPR,S_PPV,S_FM,S_OL,S_SC] = SegEva_Object(GT_PixelIdx, BW_PixelIdx)
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    S_TP = 0;
    S_FN = 0;
    
    for i = 1:length(GT_PixelIdx)
        pixelFound = false;
        for j = 1:length(BW_PixelIdx)
            if(GT_PixelIdx(i) == BW_PixelIdx(j))
                S_TP = S_TP + 1;
                pixelFound = true;
                break;
            end
        end
        if(pixelFound == false)
            S_FN = S_FN + 1;
        end
    end
    S_FP = length(BW_PixelIdx) - S_TP;
    
    % True Positive Rate or Sensitivity or Recall
    S_TPR = S_TP / ( S_TP + S_FN );
    % Positive Predictive Value or Precision or Confidence
    S_PPV = S_TP / ( S_TP + S_FP );
    % F-Measure or F-Score
    S_FM = 2 * S_TPR * S_PPV / ( S_TPR + S_PPV );
    % Overlap Ration (Jaccard Coefficient
    S_OL = S_TP / ( S_TP + S_FP + S_FN );
    % Similarity Coefficient
    S_SC = 1 - ( abs( S_FN - S_FP ) / ( 2 * S_TP + S_FN + S_FP ));
end

