function [Result] = CreateImageEvaluationTable( Files )

    for i = 1:length(Files)
        [~,imagename,~] = fileparts(Files(i).name);
        ImageName(i) = cellstr(imagename);
    end
    ImageName = ImageName';

    D_TP = zeros(length(Files),1);
    D_FN = zeros(length(Files),1);
    D_FP = zeros(length(Files),1);
    D_TPR = zeros(length(Files),1);
    D_PPV = zeros(length(Files),1);
    D_FM = zeros(length(Files),1);
    
    S_TPR = zeros(length(Files),1);
    S_PPV = zeros(length(Files),1);
    S_FM = zeros(length(Files),1);
    S_TNR = zeros(length(Files),1);
    S_FPR = zeros(length(Files),1);
    S_AUC = zeros(length(Files),1);
    S_Acc = zeros(length(Files),1);
    S_OL = zeros(length(Files),1);
    S_Dic = zeros(length(Files),1);
    S_Dic_N = zeros(length(Files),1);
    S_Kappa = zeros(length(Files),1);
    S_SC = zeros(length(Files),1);

    S_RI = zeros(length(Files),1);
    S_ARI = zeros(length(Files),1);
    S_GCE = zeros(length(Files),1);
    S_MI = zeros(length(Files),1);
    S_VI = zeros(length(Files),1);
                     
    Result = table(ImageName, D_TP, D_FN, D_FP, D_TPR, D_PPV, D_FM, ...
                    S_TPR, S_PPV, S_FM, S_TNR, S_FPR, S_AUC, S_Acc, ...
                    S_OL, S_Dic, S_Dic_N, S_Kappa, S_SC, S_RI, ...
                    S_ARI, S_GCE, S_MI, S_VI );
    
end
