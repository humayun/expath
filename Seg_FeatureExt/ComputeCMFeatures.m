function [ GLCMFeatures ] = ComputeCMFeatures( RGB, BW, NoOfChannels, GrayLevels, FeaturesPath, ImageName)
%ComputeCMFeatures compute Co-occurrence features in seven color 
%channels (Red, Green, Blue, HSV(V), Lab(L), H&E(H) and BR
%
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    if(nargin < 4)
        GrayLevels = 64;
    end
    
     myFcn=@(x)graycomatrix(x,'Offset',[1 0;0 1;1 1;-1 1],'Symmetric',...
         true, 'NumLevels', GrayLevels);

    [x,y]=size(BW);
    [L,num] = bwlabel(BW);
    bb = regionprops(L,'BoundingBox'); 
    bb = floor(cat(1,bb.BoundingBox));

    RGB = uint8(256*mat2gray(RGB));
    GS_or_RGB = size(RGB,3);
    
    if (GS_or_RGB == 1)
        GLCMFeatures = struct( ...
        'Autocorrelation_GS',0, 'CorrelationP_GS',0, 'Contrast_GS',0, ...
        'ClusterShade_GS',0, 'ClusterProminence_GS',0, 'Energy_GS',0, ...
        'Entropy_GS',0, 'HomogeneityP_GS',0, 'InverseDiffNormalized_GS',0, ...
        'InverseDiffMomentNormalized_GS',0, 'Dissimilarity_GS',0, ...
        'MaxProbability_GS',0, 'InfoMeasureCorr1_GS',0, ... 
        'InformeasureCorr2_GS',0);

        for i = 1:num
            obj=RGB;
            obj(~(L==i)) = NaN;
            top_x = bb(i,2);
            if(top_x<1)
                top_x = 1;
            end
            top_y = bb(i,1);
            if(top_y<1)
                top_y = 1;
            end
            width = top_x+bb(i,4);
            if(width > x)
                width = x;
            end
            height = top_y+bb(i,3);
            if(height > y)
                height = y;
            end
            patch = obj(top_x:width,top_y:height);
            glcm = myFcn(patch);
            glcmFeatures = GLCM_Features(glcm,0);
            GLCMFeatures.Autocorrelation_GS(i,1) = mean(glcmFeatures.autoc);
            GLCMFeatures.CorrelationP_GS(i,1) = mean(glcmFeatures.corrp);
            GLCMFeatures.Contrast_GS(i,1) = mean(glcmFeatures.contr);
            GLCMFeatures.ClusterShade_GS(i,1) = mean(glcmFeatures.cshad);
            GLCMFeatures.ClusterProminence_GS(i,1) = mean(glcmFeatures.cprom);
            GLCMFeatures.Energy_GS(i,1) = mean(glcmFeatures.energ);
            GLCMFeatures.Entropy_GS(i,1) = mean(glcmFeatures.entro);
            GLCMFeatures.HomogeneityP_GS(i,1) = mean(glcmFeatures.homop);
            GLCMFeatures.InverseDiffNormalized_GS(i,1) = mean(glcmFeatures.indnc);
            GLCMFeatures.InverseDiffMomentNormalized_GS(i,1) = mean(glcmFeatures.idmnc);
            GLCMFeatures.Dissimilarity_GS(i,1) = mean(glcmFeatures.dissi);
            GLCMFeatures.MaxProbability_GS(i,1) = mean(glcmFeatures.maxpr);
            GLCMFeatures.InfoMeasureCorr1_GS(i,1) = mean(glcmFeatures.inf1h);
            GLCMFeatures.InformeasureCorr2_GS(i,1) = mean(glcmFeatures.inf2h);
        end
        
    elseif(GS_or_RGB == 3)
        
        if(NoOfChannels == 4 && isnumeric(NoOfChannels))
            GLCMFeatures = struct( ...
            'Autocorrelation_R',0, 'CorrelationP_R',0, 'Contrast_R',0, ...
            'ClusterShade_R',0, 'ClusterProminence_R',0, 'Energy_R',0, ...
            'Entropy_R',0, 'HomogeneityP_R',0, 'InverseDiffNormalized_R',0, ...
            'InverseDiffMomentNormalized_R',0, 'Dissimilarity_R',0, ...
            'MaxProbability_R',0, 'InfoMeasureCorr1_R',0, ... 
            'InformeasureCorr2_R',0, ...
            'Autocorrelation_V',0, 'CorrelationP_V',0, 'Contrast_V',0, ...
            'ClusterShade_V',0, 'ClusterProminence_V',0, 'Energy_V',0, ...
            'Entropy_V',0, 'HomogeneityP_V',0, 'InverseDiffNormalized_V',0, ...
            'InverseDiffMomentNormalized_V',0, 'Dissimilarity_V',0, ...
            'MaxProbability_V',0, 'InfoMeasureCorr1_V',0, ... 
            'InformeasureCorr2_V',0, ...
            'Autocorrelation_L',0, 'CorrelationP_L',0, 'Contrast_L',0, ...
            'ClusterShade_L',0, 'ClusterProminence_L',0, 'Energy_L',0, ...
            'Entropy_L',0, 'HomogeneityP_L',0, 'InverseDiffNormalized_L',0, ...
            'InverseDiffMomentNormalized_L',0, 'Dissimilarity_L',0, ...
            'MaxProbability_L',0, 'InfoMeasureCorr1_L',0, ... 
            'InformeasureCorr2_L',0, ...
            'Autocorrelation_H',0, 'CorrelationP_H',0, 'Contrast_H',0, ...
            'ClusterShade_H',0, 'ClusterProminence_H',0, 'Energy_H',0, ...
            'Entropy_H',0, 'HomogeneityP_H',0, 'InverseDiffNormalized_H',0, ...
            'InverseDiffMomentNormalized_H',0, 'Dissimilarity_H',0, ...
            'MaxProbability_H',0, 'InfoMeasureCorr1_H',0, ... 
            'InformeasureCorr2_H',0);

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch);
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_R(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_R(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_R(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_R(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_R(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_R(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_R(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_R(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_R(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_R(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_R(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_R(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_R(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_R(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 2 - V (HSV) Channel
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_V(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_V(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_V(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_V(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_V(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_V(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_V(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_V(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_V(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_V(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_V(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_V(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_V(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_V(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 3 - L (Lab) Channel
            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_L(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_L(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_L(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_L(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_L(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_L(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_L(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_L(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_L(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_L(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_L(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_L(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_L(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_L(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 4 - H (H&E) Channel
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_H(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_H(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_H(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_H(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_H(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_H(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_H(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_H(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_H(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_H(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_H(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_H(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_H(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_H(i,1) = mean(glcmFeatures.inf2h);
            end

        elseif(NoOfChannels == 7 && isnumeric(NoOfChannels))

            GLCMFeatures = struct( ...
            'Autocorrelation_R',0, 'CorrelationP_R',0, 'Contrast_R',0, ...
            'ClusterShade_R',0, 'ClusterProminence_R',0, 'Energy_R',0, ...
            'Entropy_R',0, 'HomogeneityP_R',0, 'InverseDiffNormalized_R',0, ...
            'InverseDiffMomentNormalized_R',0, 'Dissimilarity_R',0, ...
            'MaxProbability_R',0, 'InfoMeasureCorr1_R',0, ... 
            'InformeasureCorr2_R',0, ...
            'Autocorrelation_G',0, 'CorrelationP_G',0, 'Contrast_G',0, ...
            'ClusterShade_G',0, 'ClusterProminence_G',0, 'Energy_G',0, ...
            'Entropy_G',0, 'HomogeneityP_G',0, 'InverseDiffNormalized_G',0, ...
            'InverseDiffMomentNormalized_G',0, 'Dissimilarity_G',0, ...
            'MaxProbability_G',0, 'InfoMeasureCorr1_G',0, ... 
            'InformeasureCorr2_G',0, ...
            'Autocorrelation_B',0, 'CorrelationP_B',0, 'Contrast_B',0, ...
            'ClusterShade_B',0, 'ClusterProminence_B',0, 'Energy_B',0, ...
            'Entropy_B',0, 'HomogeneityP_B',0, 'InverseDiffNormalized_B',0, ...
            'InverseDiffMomentNormalized_B',0, 'Dissimilarity_B',0, ...
            'MaxProbability_B',0, 'InfoMeasureCorr1_B',0, ... 
            'InformeasureCorr2_B',0, ...
            'Autocorrelation_V',0, 'CorrelationP_V',0, 'Contrast_V',0, ...
            'ClusterShade_V',0, 'ClusterProminence_V',0, 'Energy_V',0, ...
            'Entropy_V',0, 'HomogeneityP_V',0, 'InverseDiffNormalized_V',0, ...
            'InverseDiffMomentNormalized_V',0, 'Dissimilarity_V',0, ...
            'MaxProbability_V',0, 'InfoMeasureCorr1_V',0, ... 
            'InformeasureCorr2_V',0, ...
            'Autocorrelation_L',0, 'CorrelationP_L',0, 'Contrast_L',0, ...
            'ClusterShade_L',0, 'ClusterProminence_L',0, 'Energy_L',0, ...
            'Entropy_L',0, 'HomogeneityP_L',0, 'InverseDiffNormalized_L',0, ...
            'InverseDiffMomentNormalized_L',0, 'Dissimilarity_L',0, ...
            'MaxProbability_L',0, 'InfoMeasureCorr1_L',0, ... 
            'InformeasureCorr2_L',0, ...
            'Autocorrelation_H',0, 'CorrelationP_H',0, 'Contrast_H',0, ...
            'ClusterShade_H',0, 'ClusterProminence_H',0, 'Energy_H',0, ...
            'Entropy_H',0, 'HomogeneityP_H',0, 'InverseDiffNormalized_H',0, ...
            'InverseDiffMomentNormalized_H',0, 'Dissimilarity_H',0, ...
            'MaxProbability_H',0, 'InfoMeasureCorr1_H',0, ... 
            'InformeasureCorr2_H',0, ...
            'Autocorrelation_Br',0, 'CorrelationP_Br',0, 'Contrast_Br',0, ...
            'ClusterShade_Br',0, 'ClusterProminence_Br',0, 'Energy_Br',0, ...
            'Entropy_Br',0, 'HomogeneityP_Br',0, 'InverseDiffNormalized_Br',0, ...
            'InverseDiffMomentNormalized_Br',0, 'Dissimilarity_Br',0, ...
            'MaxProbability_Br',0, 'InfoMeasureCorr1_Br',0, ... 
            'InformeasureCorr2_Br',0 );

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch);
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_R(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_R(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_R(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_R(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_R(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_R(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_R(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_R(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_R(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_R(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_R(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_R(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_R(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_R(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 2 - Green Channel
            channel = RGB(:,:,2);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_G(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_G(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_G(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_G(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_G(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_G(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_G(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_G(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_G(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_G(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_G(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_G(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_G(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_G(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 3 - Blue Channel
            channel = RGB(:,:,3);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_B(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_B(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_B(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_B(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_B(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_B(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_B(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_B(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_B(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_B(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_B(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_B(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_B(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_B(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 4 - V (HSV) Channel
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_V(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_V(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_V(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_V(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_V(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_V(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_V(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_V(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_V(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_V(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_V(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_V(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_V(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_V(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 5 - L (Lab) Channel
            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_L(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_L(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_L(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_L(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_L(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_L(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_L(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_L(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_L(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_L(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_L(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_L(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_L(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_L(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 6 - H (H&E) Channel
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_H(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_H(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_H(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_H(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_H(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_H(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_H(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_H(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_H(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_H(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_H(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_H(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_H(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_H(i,1) = mean(glcmFeatures.inf2h);
            end

            %% 7 - BlueRatio Image
            channel= uint8(GrayLevels*mat2gray(ComputeBlueRatio(RGB)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_Br(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_Br(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_Br(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_Br(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_Br(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_Br(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_Br(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_Br(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_Br(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_Br(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_Br(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_Br(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_Br(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_Br(i,1) = mean(glcmFeatures.inf2h);
            end
            
        elseif(NoOfChannels == 'R' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_R',0, 'CorrelationP_R',0, 'Contrast_R',0, ...
            'ClusterShade_R',0, 'ClusterProminence_R',0, 'Energy_R',0, ...
            'Entropy_R',0, 'HomogeneityP_R',0, 'InverseDiffNormalized_R',0, ...
            'InverseDiffMomentNormalized_R',0, 'Dissimilarity_R',0, ...
            'MaxProbability_R',0, 'InfoMeasureCorr1_R',0, ... 
            'InformeasureCorr2_R',0);

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch);
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_R(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_R(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_R(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_R(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_R(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_R(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_R(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_R(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_R(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_R(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_R(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_R(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_R(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_R(i,1) = mean(glcmFeatures.inf2h);
            end

        elseif(NoOfChannels == 'G' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_G',0, 'CorrelationP_G',0, 'Contrast_G',0, ...
            'ClusterShade_G',0, 'ClusterProminence_G',0, 'Energy_G',0, ...
            'Entropy_G',0, 'HomogeneityP_G',0, 'InverseDiffNormalized_G',0, ...
            'InverseDiffMomentNormalized_G',0, 'Dissimilarity_G',0, ...
            'MaxProbability_G',0, 'InfoMeasureCorr1_G',0, ... 
            'InformeasureCorr2_G',0);

            channel = RGB(:,:,2);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_G(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_G(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_G(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_G(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_G(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_G(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_G(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_G(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_G(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_G(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_G(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_G(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_G(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_G(i,1) = mean(glcmFeatures.inf2h);
            end

        elseif(NoOfChannels == 'B' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_B',0, 'CorrelationP_B',0, 'Contrast_B',0, ...
            'ClusterShade_B',0, 'ClusterProminence_B',0, 'Energy_B',0, ...
            'Entropy_B',0, 'HomogeneityP_B',0, 'InverseDiffNormalized_B',0, ...
            'InverseDiffMomentNormalized_B',0, 'Dissimilarity_B',0, ...
            'MaxProbability_B',0, 'InfoMeasureCorr1_B',0, ... 
            'InformeasureCorr2_B',0);
        
            channel = RGB(:,:,3);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_B(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_B(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_B(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_B(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_B(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_B(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_B(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_B(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_B(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_B(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_B(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_B(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_B(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_B(i,1) = mean(glcmFeatures.inf2h);
            end
            
        elseif(NoOfChannels == 'V' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_V',0, 'CorrelationP_V',0, 'Contrast_V',0, ...
            'ClusterShade_V',0, 'ClusterProminence_V',0, 'Energy_V',0, ...
            'Entropy_V',0, 'HomogeneityP_V',0, 'InverseDiffNormalized_V',0, ...
            'InverseDiffMomentNormalized_V',0, 'Dissimilarity_V',0, ...
            'MaxProbability_V',0, 'InfoMeasureCorr1_V',0, ... 
            'InformeasureCorr2_V',0);
        
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_V(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_V(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_V(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_V(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_V(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_V(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_V(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_V(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_V(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_V(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_V(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_V(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_V(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_V(i,1) = mean(glcmFeatures.inf2h);
            end
            
        elseif(NoOfChannels == 'L' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_L',0, 'CorrelationP_L',0, 'Contrast_L',0, ...
            'ClusterShade_L',0, 'ClusterProminence_L',0, 'Energy_L',0, ...
            'Entropy_L',0, 'HomogeneityP_L',0, 'InverseDiffNormalized_L',0, ...
            'InverseDiffMomentNormalized_L',0, 'Dissimilarity_L',0, ...
            'MaxProbability_L',0, 'InfoMeasureCorr1_L',0, ... 
            'InformeasureCorr2_L',0);

            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_L(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_L(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_L(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_L(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_L(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_L(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_L(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_L(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_L(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_L(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_L(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_L(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_L(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_L(i,1) = mean(glcmFeatures.inf2h);
            end
            
        elseif(NoOfChannels == 'H' && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_H',0, 'CorrelationP_H',0, 'Contrast_H',0, ...
            'ClusterShade_H',0, 'ClusterProminence_H',0, 'Energy_H',0, ...
            'Entropy_H',0, 'HomogeneityP_H',0, 'InverseDiffNormalized_H',0, ...
            'InverseDiffMomentNormalized_H',0, 'Dissimilarity_H',0, ...
            'MaxProbability_H',0, 'InfoMeasureCorr1_H',0, ... 
            'InformeasureCorr2_H',0);

            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_H(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_H(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_H(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_H(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_H(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_H(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_H(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_H(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_H(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_H(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_H(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_H(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_H(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_H(i,1) = mean(glcmFeatures.inf2h);
            end
            
        elseif(strcmp(NoOfChannels,'BR') && ~isnumeric(NoOfChannels))
            
            GLCMFeatures = struct( ...
            'Autocorrelation_Br',0, 'CorrelationP_Br',0, 'Contrast_Br',0, ...
            'ClusterShade_Br',0, 'ClusterProminence_Br',0, 'Energy_Br',0, ...
            'Entropy_Br',0, 'HomogeneityP_Br',0, 'InverseDiffNormalized_Br',0, ...
            'InverseDiffMomentNormalized_Br',0, 'Dissimilarity_Br',0, ...
            'MaxProbability_Br',0, 'InfoMeasureCorr1_Br',0, ... 
            'InformeasureCorr2_Br',0 );

            channel= uint8(GrayLevels*mat2gray(ComputeBlueRatio(RGB)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                glcm = myFcn(patch); 
                glcmFeatures = GLCM_Features(glcm,0);
                GLCMFeatures.Autocorrelation_Br(i,1) = mean(glcmFeatures.autoc);
                GLCMFeatures.CorrelationP_Br(i,1) = mean(glcmFeatures.corrp);
                GLCMFeatures.Contrast_Br(i,1) = mean(glcmFeatures.contr);
                GLCMFeatures.ClusterShade_Br(i,1) = mean(glcmFeatures.cshad);
                GLCMFeatures.ClusterProminence_Br(i,1) = mean(glcmFeatures.cprom);
                GLCMFeatures.Energy_Br(i,1) = mean(glcmFeatures.energ);
                GLCMFeatures.Entropy_Br(i,1) = mean(glcmFeatures.entro);
                GLCMFeatures.HomogeneityP_Br(i,1) = mean(glcmFeatures.homop);
                GLCMFeatures.InverseDiffNormalized_Br(i,1) = mean(glcmFeatures.indnc);
                GLCMFeatures.InverseDiffMomentNormalized_Br(i,1) = mean(glcmFeatures.idmnc);
                GLCMFeatures.Dissimilarity_Br(i,1) = mean(glcmFeatures.dissi);
                GLCMFeatures.MaxProbability_Br(i,1) = mean(glcmFeatures.maxpr);
                GLCMFeatures.InfoMeasureCorr1_Br(i,1) = mean(glcmFeatures.inf1h);
                GLCMFeatures.InformeasureCorr2_Br(i,1) = mean(glcmFeatures.inf2h);
            end
        end
                                
    end
    
     if(nargin == 6)
         writetable(struct2table(GLCMFeatures, strcat(FeaturesPath, ImageName, '_CMFeatures.csv')));
     end
     GLCMFeatures = struct2table(GLCMFeatures);
end
