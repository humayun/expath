function[S_TPR,S_PPV,S_FM,S_TNR,S_FPR,S_AUC,S_Acc,S_OL,S_Dic,S_Kappa,S_SC, ...
    S_RI,S_ARI,S_GCE,S_MI,S_VI] = SegEva_Image(GT,BW)

% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    [TP,FN,FP,TN,n_GT,n_BW,A,B,C,D] = ComputeContingencyTable(GT,BW);
    
    % True Positive Rate or Sensitivity or Recall
    S_TPR = TP / ( TP + FN );
    % Positive Predictive Value or Precision or Confidence
    S_PPV = TP / ( TP + FP );
    % F-Measure or F-Score
    S_FM = 2 * S_TPR * S_PPV / ( S_TPR + S_PPV );
    % True Negative Rate or Specificity
    S_TNR = TN / ( TN + FP );
    % False Positive Rate or 
    S_FPR = 1 - S_TNR;
    % Area under the Curve of Recieve Operating Curve (ROC)
    S_AUC = ( S_TPR - S_FPR + 1 ) / 2;
    % Accuracy
    S_Acc = ( TP + TN ) / ( TP + TN + FP + FN ); 
    % Overlap Ration (Jaccard Coefficient)
    S_OL = TP / ( TP + FN + FP );
    % Dice Index
    S_Dic = 2 * TP / ( n_GT + n_BW );

    % Cohin Kappa
    agreement = TP + TN;
    chance_0 = ( TN + FN ) * ( TN + FP );
    chance_1 = ( FP + TP ) * ( FN + TP );
    N = ( TP + TN + FP + FN);
    chance = (chance_0+chance_1) / N;
    S_Kappa = (agreement - chance) / (N - chance);
    
    % Similarity Coefficient
    S_SC = 1 - (abs(FN-FP) / (2*TP + FN + FP ));
    
    % Rand Index
    S_RI = ( A + D ) / ( A + B + C + D );
        
    % Adjusted Rand Index
    X1 = A - ( ( A+C)*(A+B) / (A+B+C+D) );
    X2 = ( (A+C) + (A+B) ) / 2;
    X3 = ( (A+C) * (A+B) ) / (A+B+C+D);
    if(X2==X3)
        S_ARI = 0;
    else
        S_ARI = X1/(X2-X3);
    end
                
    % Global Consistency Error
    e1 = ( FN * ( FN + 2*TP ) / ( TP + FN ) + FP * ( FP + 2 * TN ) / ( TN + FP ) ) / N;
    e2 = ( FP * ( FP + 2*TP ) / ( TP + FP ) + FN * ( FN + 2 * TN ) / ( TN + FN ) ) / N;
    S_GCE = min(e1,e2);
    
    % Mutual Information
    row1 = TN + FN;
    row2 = FP + TP;
    H1 = - ( ( row1 / N ) * log2( row1 / N ) + ( row2 / N ) * log2( row2 / N ) );
    col1 = TN + FP;
    col2 = FN + TP;
    H2 = - ( ( col1 / N ) * log2( col1 / N ) + ( col2 / N ) * log2( col2 / N ) );
    if( TN == 0)
        p00 = 1;
    else
        p00 = TN / N;
    end
    if( FN == 0)
        p01 = 1;
    else
        p01 = FN / N;
    end
    if( FP == 0)
        p10 = 1;
    else
        p10 = FP / N;
    end
    if( TP == 0)
        p11 = 1;
    else
        p11 = TP / N;
    end
    H12 = - ( ( TN / N ) * log2( p00 ) + ( FN / N ) * log2( p01 ) + ...
            ( FP / N ) * log2( p10 ) + ( TP / N ) * log2( p11 ) );
    S_MI = H1 + H2 - H12;
    
    % Variation of Information
    S_VI = H1 + H2 - (2 * S_MI);
    
end
