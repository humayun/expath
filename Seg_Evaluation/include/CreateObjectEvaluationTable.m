function [Result] = CreateObjectEvaluationTable( NoOfObjects )

    N_Found = zeros(NoOfObjects,1);
    S_TP = zeros(NoOfObjects,1);
    S_FN = zeros(NoOfObjects,1);
    S_FP = zeros(NoOfObjects,1);
    S_TPR = zeros(NoOfObjects,1);
    S_PPV = zeros(NoOfObjects,1);
    S_FM = zeros(NoOfObjects,1);
    S_OL = zeros(NoOfObjects,1);
    S_SC = zeros(NoOfObjects,1);
    S_Dic = zeros(NoOfObjects,1);
                     
    Result = table(N_Found,S_TP, S_FN, S_FP, S_TPR, S_PPV, S_FM, S_OL, S_SC, S_Dic);
    
end