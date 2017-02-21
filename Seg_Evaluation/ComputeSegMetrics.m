% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

% addpath('Include');
% addpath('LabelMeToolbox');
% addpath('LabelMeToolbox/3Dtools');
% addpath('LabelMeToolbox/ADMINtools');
% addpath('LabelMeToolbox/SUNdatabase');
% addpath('LabelMeToolbox/XMLtools');
% addpath('LabelMeToolbox/compatibility');
% addpath('LabelMeToolbox/features');
% addpath('LabelMeToolbox/imagemanipulation');
% addpath('LabelMeToolbox/install');
% addpath('LabelMeToolbox/main');
% addpath('LabelMeToolbox/objectdetection');
% addpath('LabelMeToolbox/old');
% addpath('LabelMeToolbox/parts');
% addpath('LabelMeToolbox/primalSVM');
% addpath('LabelMeToolbox/querytools');
% addpath('LabelMeToolbox/utils');
% addpath('LabelMeToolbox/video');
% addpath('LabelMeToolbox/wordnet');

%% 
%*** This program required LabelMe Toolbox, This program download data
%either from LabelMe website or open annotation images from already
%downloaded dataset.

%% Downloadiing User Database (Images + Annotations) from LabelMe Website
% Define the root folder for the images
HOMEIMAGES = 'Expanded/';
HOMEANNOTATIONS = strcat(HOMEIMAGES,'Annotation'); 
D = LMdatabase(HOMEANNOTATIONS);    % Reading the index
HOMEMASKS = strcat(HOMEIMAGES,'Masks/'); 
mkdir(HOMEMASKS);
% Database name with path on LabelMe website
%DatabaseName ={'users/hirshad/becklab_perexpanded'};
%DatabaseName ={'users/hirshad/becklab_expanded'};
%LMinstall(DatabaseName, HOMEIMAGES,HOMEANNOTATIONS);

%% Reading Automated Method Segmentation files
BWPath = strcat(HOMEIMAGES,'AutomatedSegmentation/');
BWExt = '.png';
Files = dir(strcat(BWPath,'*',BWExt));

%% Create empty Table for storing metrics
IR = CreateImageEvaluationTable(Files);

%% Extracting Mask for all objects
for i=1:length(D)
    % Reading Ground Truth from LabelMe Dataset
    [mask, ~] = LMobjectmask(D(i).annotation,HOMEIMAGES);
    GT = logical(colorSegments(mask));
    GT = max(max(GT(:,:,1), GT(:,:,2)), GT(:,:,3));
    GT_F = regionprops(logical(GT),'Area','Centroid','BoundingBox',...
                                            'PixelIdxList','PixelList');
    [~,~,N] = size(mask);
    NR = CreateObjectEvaluationTable(N);
    
    % Reading Automated Segmentation Image
    [~,ImageName,~] = fileparts(D(i).annotation.filename);
    BW = imread( strcat( BWPath, ImageName, '_binary', BWExt));
    BW_F = regionprops(logical(BW),'Area','Centroid','BoundingBox',...
                                            'PixelIdxList','PixelList');
    % A array of (intially) zeros for correcting detecting segmented nuclei 
    BW_NucleiFound = zeros(length(BW_F),1);
    
    % loop for each nuclei in ground truth mask
    for j = 1:N
        GT_SN_F = regionprops(logical(mask(:,:,j)),'Area','Centroid',...
                                'BoundingBox','PixelIdxList','PixelList');
        if(~isempty(GT_SN_F))
            found = false;
            % loop for each nuclei in segmentet mask
            for k=1:length(BW_F)
                % if the centroids of segmented nuclei in all pixels of GT nuclei
                if(sum(GT_SN_F(1).PixelList(:,1) == int16(BW_F(k).Centroid(1))) > 0 && ...
                       sum(GT_SN_F(1).PixelList(:,2) == int16(BW_F(k).Centroid(2))) > 0 && ...
                       BW_NucleiFound(k) == 0)
                   
                    BW_NucleiFound(k) = 1;
                    found = true;
                    NR.N_Found(j) = 1;
                    IR.D_TP(i) = IR.D_TP(i) + 1;

                    % Segmentation Metrics at Nuclei Level
                    [NR.S_TP(j),NR.S_FN(j),NR.S_FP(j),NR.S_TPR(j),NR.S_PPV(j), ...
                        NR.S_FM(j),NR.S_OL(j),NR.S_SC(j)] ... 
                        = SegEva_Object( GT_SN_F(1).PixelIdxList, BW_F(k).PixelIdxList);

                    NR.S_Dic(j) = 2 * NR.S_TP(j) / ( GT_SN_F(1).Area + BW_F(k).Area );

                    % Segmentation Metric at Image Level
                    alpha = GT_SN_F(1).Area / sum(cat(1,GT_F.Area));
                    beta = BW_F(k).Area / sum(cat(1,BW_F.Area));
                    IR.S_Dic_N(i) = IR.S_Dic_N(i) + ...
                            ( ( (alpha * NR.S_Dic(j)) + (beta * NR.S_Dic(j)) ) / 2 );

                    break;
                end
            end
            if(~found)
                IR.D_FN(i) = IR.D_FN(i) + 1;
                NR.N_Found(j) = -1;
            end
        end
    end
    
    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( BWPath, ImageName, '_NucleiMetrics.csv'));
    
    % Detection Metrics at Image Level
    IR.D_FP(i) = length(BW_F) - IR.D_TP(i);
    IR.D_TPR(i) = IR.D_TP(i) / (IR.D_TP(i) + IR.D_FN(i));
    IR.D_PPV(i) = IR.D_TP(i) / (IR.D_TP(i) + IR.D_FP(i));
    IR.D_FM(i) = 2 * IR.D_TPR(i) * IR.D_PPV(i) / ...
                        (IR.D_TPR(i) + IR.D_PPV(i));

    % Segmentation Metrics at Image Level
    [ IR.S_TPR(i),IR.S_PPV(i),IR.S_FM(i),IR.S_TNR(i),IR.S_FPR(i), ...
        IR.S_AUC(i),IR.S_Acc(i),IR.S_OL(i),IR.S_Dic(i),IR.S_Kappa(i), ...
        IR.S_SC(i),IR.S_RI(i),IR.S_ARI(i),IR.S_GCE(i), ...
        IR.S_MI(i),IR.S_VI(i) ] = SegEva_Image(GT,BW);
end

% Saving each image results (Image Level Metrics)
writetable(IR,strcat( BWPath, 'ImageDetectionAndSegmentationMetrics.csv'));
