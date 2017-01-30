function [ IntensityFeatures ] = ComputeIntensityFeatures( RGB, BW, NoOfChannels, FeaturesPath, ImageName)
% ComputeIntensityFeatures compute intesnity features (Mean, Median, MAD,
% Standard Deviation, Inter-Quantile Range, Skewness and Kurtosis in 
% grayscale or RGB image. In case of RGB color image, it compute features
% in selected color channels (Red, Green, Blue, HSV(V), Lab(L), H&E(H) and BR
% If NoOfChannels is number then either 4 or 7, then this function compute 
% features in 4 selected or 7 selected channels. If NoOfChannels is
% character, then it compute features only in that channels.
%
% author: Humayun Irshad (12/2015)

%%
    nucleiStats = regionprops( logical(BW), 'PixelIdxList' );

    GS_or_RGB = size(RGB,3);
    if (GS_or_RGB == 1)
        IntensityFeatures = struct( ...
          'Min_GS',0,'Max_GS',0,'Mean_GS',0,'Median_GS',0,'MAD_GS',0,'SD_GS',0,'IQR_GS',0, ...
          'Skewness_GS',0,'Kurtosis_GS',0);

        for i = 1:numel(nucleiStats)
            nucleiPixelIntensities = RGB( nucleiStats(i).PixelIdxList );
            IntensityFeatures.Min_GS(i,1) = min( nucleiPixelIntensities );
            IntensityFeatures.Max_GS(i,1) = max( nucleiPixelIntensities );
            IntensityFeatures.Mean_GS(i,1) = mean(nucleiPixelIntensities);
            IntensityFeatures.Median_GS(i,1) = median(nucleiPixelIntensities);
            IntensityFeatures.MAD_GS(i,1) = mad(double(nucleiPixelIntensities));
            IntensityFeatures.SD_GS(i,1) = std(double(nucleiPixelIntensities));
            IntensityFeatures.IQR_GS(i,1) = iqr(double(nucleiPixelIntensities));
            IntensityFeatures.Skewness_GS(i,1) = skewness(double(nucleiPixelIntensities));
            %IntensityFeatures.SkewnessRH_GS(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
            IntensityFeatures.Kurtosis_GS(i,1) = kurtosis(double(nucleiPixelIntensities));
            %IntensityFeatures.KurtosisRC_GS(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
        end
        
    elseif(GS_or_RGB == 3)
        if(NoOfChannels == 4 && isnumeric(NoOfChannels))
            IntensityFeatures = struct( ...
               'Min_R',0,'Max_R',0,'Mean_R',0,'Median_R',0,'MAD_R',0,'SD_R',0,'IQR_R',0, ...
              'Skewness_R',0,'Kurtosis_R',0, ...
              'Min_V',0,'Max_V',0,'Mean_V',0,'Median_V',0,'MAD_V',0,'SD_V',0,'IQR_V',0, ...
              'Skewness_V',0,'Kurtosis_V',0, ...
              'Min_L',0,'Max_L',0,'Mean_L',0,'Median_L',0,'MAD_L',0,'SD_L',0,'IQR_L',0, ...
              'Skewness_L',0,'Kurtosis_L',0, ...
              'Min_H',0,'Max_H',0,'Mean_H',0,'Median_H',0,'MAD_H',0,'SD_H',0,'IQR_H',0, ...
              'Skewness_H',0,'Kurtosis_H',0);

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_R(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_R(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_R(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_R(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_R(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_R(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_R(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_R(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_R(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_R(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_R(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 2 - V (HSV) Color Channel
            channel = rgb2hsv(RGB);
            channel = uint8(255*mat2gray(channel(:,:,3)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_V(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_V(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_V(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_V(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_V(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_V(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_V(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_V(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_V(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_V(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_V(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 3 - L (Lab) Color Channel
            channel = rgb2lab(RGB);
            channel = uint8(255*mat2gray(channel(:,:,1)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_L(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_L(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_L(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_L(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_L(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_L(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_L(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_L(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_L(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_L(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_L(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 4 - H (H&E) Color Deconvolution
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(256*mat2gray(channel));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_H(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_H(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_H(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_H(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_H(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_H(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_H(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_H(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_H(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_H(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_H(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 7 && isnumeric(NoOfChannels))
            IntensityFeatures = struct( ...
              'Min_R',0,'Max_R',0,'Mean_R',0,'Median_R',0,'MAD_R',0,'SD_R',0,'IQR_R',0, ...
              'Skewness_R',0,'Kurtosis_R',0, ...
              'Min_G',0,'Max_G',0,'Mean_G',0,'Median_G',0,'MAD_G',0,'SD_G',0,'IQR_G',0, ...
              'Skewness_G',0,'Kurtosis_G',0, ...
              'Min_B',0,'Max_B',0,'Mean_B',0,'Median_B',0,'MAD_B',0,'SD_B',0,'IQR_B',0, ...
              'Skewness_B',0,'Kurtosis_B',0, ...
              'Min_V',0,'Max_V',0,'Mean_V',0,'Median_V',0,'MAD_V',0,'SD_V',0,'IQR_V',0, ...
              'Skewness_V',0,'Kurtosis_V',0, ...
              'Min_L',0,'Max_L',0,'Mean_L',0,'Median_L',0,'MAD_L',0,'SD_L',0,'IQR_L',0, ...
              'Skewness_L',0,'Kurtosis_L',0, ...
              'Min_H',0,'Max_H',0,'Mean_H',0,'Median_H',0,'MAD_H',0,'SD_H',0,'IQR_H',0, ...
              'Skewness_H',0,'Kurtosis_H',0, ...
              'Min_Br',0,'Max_Br',0,'Mean_Br',0,'Median_Br',0,'MAD_Br',0,'SD_Br',0,'IQR_Br',0, ...
              'Skewness_Br',0,'Kurtosis_Br',0);

            for i = 1:numel(nucleiStats)
            %% 1 - Red Channel
                channel = RGB(:,:,1);
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_R(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_R(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_R(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_R(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_R(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_R(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_R(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_R(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_R(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_R(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_R(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));

            %% 2 - Greeen Channel
                channel = RGB(:,:,2);
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_G(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_G(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_G(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_G(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_G(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_G(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_G(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_G(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_G(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_G(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_G(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));

            %% 3 - Blue Channel
                channel = RGB(:,:,3);
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_B(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_B(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_B(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_B(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_B(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_B(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_B(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_B(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_B(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_B(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_B(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 4 - V (HSV) Color Channel
            channel = rgb2hsv(RGB);
            channel = uint8(255*mat2gray(channel(:,:,3)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_V(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_V(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_V(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_V(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_V(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_V(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_V(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_V(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_V(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_V(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_V(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 5 - L (Lab) Color Channel
            channel = rgb2lab(RGB);
            channel = uint8(255*mat2gray(channel(:,:,1)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_L(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_L(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_L(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_L(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_L(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_L(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_L(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_L(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_L(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_L(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_L(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 6 - H (H&E) Color Deconvolution
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(256*mat2gray(channel));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_H(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_H(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_H(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_H(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_H(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_H(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_H(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_H(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_H(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_H(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_H(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end

            %% 7 - BlueRatio Image Channel
            channel= ComputeBlueRatio(RGB);
            channel = uint8(256*mat2gray(channel));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_Br(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_Br(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_Br(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_Br(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_Br(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_Br(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_Br(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_Br(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_Br(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_Br(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_Br(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'R' && ~isnumeric(NoOfChannels))
            %% Red Channel
            IntensityFeatures = struct( ...
              'Min_R',0,'Max_R',0,'Mean_R',0,'Median_R',0,'MAD_R',0,'SD_R',0,'IQR_R',0, ...
              'Skewness_R',0,'Kurtosis_R',0);
            
            channel = RGB(:,:,1);
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_R(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_R(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_R(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_R(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_R(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_R(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_R(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_R(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_R(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_R(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_R(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'G' && ~isnumeric(NoOfChannels))
            %% Green Channel
            IntensityFeatures = struct( ...
              'Min_G',0,'Max_G',0,'Mean_G',0,'Median_G',0,'MAD_G',0,'SD_G',0,'IQR_G',0, ...
              'Skewness_G',0,'Kurtosis_G',0);

            channel = RGB(:,:,2);
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_G(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_G(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_G(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_G(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_G(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_G(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_G(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_G(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_G(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_G(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_G(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'B' && ~isnumeric(NoOfChannels))
            %% Blue Channel
            IntensityFeatures = struct( ...
              'Min_B',0,'Max_B',0,'Mean_B',0,'Median_B',0,'MAD_B',0,'SD_B',0,'IQR_B',0, ...
              'Skewness_B',0,'Kurtosis_B',0);

            channel = RGB(:,:,3);
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_B(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_B(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_B(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_B(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_B(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_B(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_B(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_B(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_B(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_B(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_B(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'V' && ~isnumeric(NoOfChannels))
            %% V (HSV) Channel
            IntensityFeatures = struct( ...
              'Min_V',0,'Max_V',0,'Mean_V',0,'Median_V',0,'MAD_V',0,'SD_V',0,'IQR_V',0, ...
              'Skewness_V',0,'Kurtosis_V',0);

            channel = rgb2hsv(RGB);
            channel = uint8(255*mat2gray(channel(:,:,3)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_V(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_V(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_V(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_V(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_V(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_V(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_V(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_V(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_V(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_V(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_V(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'L' && ~isnumeric(NoOfChannels))
            %% L (Lab) Channel
            IntensityFeatures = struct( ...
              'Min_L',0,'Max_L',0,'Mean_L',0,'Median_L',0,'MAD_L',0,'SD_L',0,'IQR_L',0, ...
              'Skewness_L',0,'Kurtosis_L',0);

            channel = rgb2lab(RGB);
            channel = uint8(255*mat2gray(channel(:,:,1)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_L(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_L(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_L(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_L(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_L(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_L(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_L(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_L(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_L(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_L(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_L(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(NoOfChannels == 'H' && ~isnumeric(NoOfChannels))
            %% H (H&E) Channel
            IntensityFeatures = struct( ...
              'Min_H',0,'Max_H',0,'Mean_H',0,'Median_H',0,'MAD_H',0,'SD_H',0,'IQR_H',0, ...
              'Skewness_H',0,'Kurtosis_H',0);
          
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
             channel = uint8(256*mat2gray(channel));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_H(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_H(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_H(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_H(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_H(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_H(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_H(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_H(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_H(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_H(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_H(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
            
        elseif(strcmp(NoOfChannels,'BR') && ~isnumeric(NoOfChannels))
            %% BR (BlueRatio) Channel
            IntensityFeatures = struct( ...
              'Min_Br',0,'Max_Br',0,'Mean_Br',0,'Median_Br',0,'MAD_Br',0,'SD_Br',0,'IQR_Br',0, ...
              'Skewness_Br',0,'Kurtosis_Br',0);

            channel= uint8(256*mat2gray(ComputeBlueRatio(RGB)));
            for i = 1:numel(nucleiStats)
                nucleiPixelIntensities = channel( nucleiStats(i).PixelIdxList );
                IntensityFeatures.Min_Br(i,1) = min(nucleiPixelIntensities);
                IntensityFeatures.Max_Br(i,1) = max(nucleiPixelIntensities);
                IntensityFeatures.Mean_Br(i,1) = mean(nucleiPixelIntensities);
                IntensityFeatures.Median_Br(i,1) = median(nucleiPixelIntensities);
                IntensityFeatures.MAD_Br(i,1) = mad(double(nucleiPixelIntensities));
                IntensityFeatures.SD_Br(i,1) = std(double(nucleiPixelIntensities));
                IntensityFeatures.IQR_Br(i,1) = iqr(double(nucleiPixelIntensities));
                IntensityFeatures.Skewness_Br(i,1) = skewness(double(nucleiPixelIntensities));
                %IntensityFeatures.SkewnessRH_Br(i,1) = skewnessRobustHinkley(double(nucleiPixelIntensities));
                IntensityFeatures.Kurtosis_Br(i,1) = kurtosis(double(nucleiPixelIntensities));
                %IntensityFeatures.KurtosisRC_Br(i,1) = kurtosisRobustCrow(double(nucleiPixelIntensities));
            end
        end
    end
    
    %% Writing features
    if(nargin == 5)
        struct2csv(IntensityFeatures,strcat(FeaturesPath,ImageName,'_IntensityFeatures.csv'));
    end        
    IntensityFeatures = struct2table(IntensityFeatures);
end

