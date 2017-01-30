function [ SumFeatures ] = SummarizeFeatures_2( Features )
%SummarizeFeatures Summarize all the features 
%   SummarizeFeatures summarize all the features, like mean, median, STD,
%   Skewness, Kurtosis, MAD, IQR, 
    
    SumFeatures = struct;
    F = table2array(Features);
    FeaturesName = fieldnames(Features);
    FeaturesName = FeaturesName(1:length(FeaturesName)-1);
    for i=1:numel(FeaturesName)
        Name = strcat('Mean_',FeaturesName{i});
        tmp = mean(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('STD_',FeaturesName{i});
        tmp = std(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    SumFeatures = struct2table(SumFeatures);
end

