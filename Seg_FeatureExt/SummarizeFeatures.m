function [ SumFeatures ] = SummarizeFeatures( Features )
%SummarizeFeatures Summarize all the features 
%   SummarizeFeatures summarize all the features, like mean, median, STD,
%   Skewness, Kurtosis, MAD, IQR, 
%
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    
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
        Name = strcat('Median_',FeaturesName{i});
        tmp = median(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('MAD_',FeaturesName{i});
        tmp = mad(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('STD_',FeaturesName{i});
        tmp = std(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('IQR_',FeaturesName{i});
        tmp = iqr(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('Skewness_',FeaturesName{i});
        tmp = skewness(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
    for i=1:numel(FeaturesName)
        Name = strcat('Kurtosis_',FeaturesName{i});
        tmp = kurtosis(F(:,i));
        SumFeatures = setfield(SumFeatures, Name, tmp);
    end
    
%     for i=1:numel(FeaturesName)
%         Name = strcat('SkewnessB_',FeaturesName{i});
%         tmp = skewnessRobustHinkley(F(:,i));
%         SumFeatures = setfield(SumFeatures, Name, tmp);
%     end
%     
%     for i=1:numel(FeaturesName)
%         Name = strcat('KurtosisB_',FeaturesName{i});
%         tmp = kurtosisRobustCrow(F(:,i));
%         SumFeatures = setfield(SumFeatures, Name, tmp);
%     end
    SumFeatures = struct2table(SumFeatures);
end

