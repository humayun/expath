% addpath('/home/hi41/WS/MatCode/Lib');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/3Dtools');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/ADMINtools');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/SUNdatabase');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/XMLtools');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/compatibility');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/features');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/imagemanipulation');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/install');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/main');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/objectdetection');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/old');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/parts');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/primalSVM');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/querytools');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/utils');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/video');
% addpath('/home/hi41/WS/MatCode/LabelMeToolbox/wordnet');

%%
%*** This program required LabelMe Toolbox, This program download data
%either from LabelMe website or open annotation images from already
%downloaded dataset.

%% Database Path to Local Drive
%HOMEIMAGES = '/groups/becklab/Expansion/PreExpanded/';
HOMEIMAGES = 'H:/Datasets/LSM/LabelMe/PreExpanded/';
HOMEANNOTATIONS = strcat(HOMEIMAGES,'Annotation'); 
D = LMdatabase(HOMEANNOTATIONS);

BWPath = strcat(HOMEIMAGES,'AutomatedSegmentation/');
BWExt = '.png';

ResultPath = strcat(HOMEIMAGES,'Results/');
mkdir(ResultPath);

OverlayPath = strcat(HOMEIMAGES,'Overlay/');
mkdir(OverlayPath);

Files = dir(strcat(HOMEIMAGES,'*.jpg'));
IR = CreateImageEvaluationTable(Files);
IR_WN = CreateImageEvaluationTable(Files);
IR_SN = CreateImageEvaluationTable(Files);

for i=1:length(D)
    % Reading Ground Truth from LabelMe Dataset
    % ^^^ LMobjectmask give error because, it try to return a 3D images with 
    % 3rd dimension equal to number of labels. For expansion data, as image
    % is very big and there are many nuclei, then due to size constraint, a
    % modification in the function to address this problem.
    [GT, ClassLabels] = LMobjectmask_2(D(i).annotation,HOMEIMAGES);
   
%% Case 1: All Nuclei
    GT_F = regionprops(GT,'Area','Centroid','BoundingBox',...
                                            'PixelIdxList','PixelList');
    N = length(GT_F);
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
        found = false;
        % loop for each nuclei in segmentet mask
        for k=1:length(BW_F)
            % if the centroids of segmented nuclei inside GT nuclei
            if(sum(GT_F(j).PixelList(:,1) == int16(BW_F(k).Centroid(1))) > 0 && ...
                   sum(GT_F(j).PixelList(:,2) == int16(BW_F(k).Centroid(2))) > 0 && ...
                   BW_NucleiFound(k) == 0)
               
                BW_NucleiFound(k) = 1;
                found = true;
                NR.N_Found(j) = 1;
                IR.D_TP(i) = IR.D_TP(i) + 1;

                % Segmentation Metrics at Nuclei Level
                [NR.S_TP(j),NR.S_FN(j),NR.S_FP(j),NR.S_TPR(j),NR.S_PPV(j), ...
                    NR.S_FM(j),NR.S_OL(j),NR.S_SC(j)] ... 
                    = SegEva_Object( GT_F(j).PixelIdxList, BW_F(k).PixelIdxList);

                NR.S_Dic(j) = 2 * NR.S_TP(j) / ( GT_F(j).Area + BW_F(k).Area );

                % Segmentation Metric at Image Level
                alpha = GT_F(j).Area / sum(cat(1,GT_F.Area));
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

    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( ResultPath, ImageName, '_Metrics_AllNuclei.csv'));
    imwrite(GT, strcat(OverlayPath,ImageName,'_AllNuclei.png'));
    
    im = imread(strcat(HOMEIMAGES,ImageName,'.jpg'));
    
    % BW 
    im_BW = imoverlay_v2(im,BW,'color',[0 1 0],'alpha',0.75,'bright');
    imwrite( im_BW, strcat(OverlayPath,ImageName,'_BW_Overlay.png') );
    im_P = imdilate(boundarymask(BW), ones(3,3));
    imwrite(imoverlay(im,im_P,'green'),strcat(OverlayPath,ImageName,'_BW_Perimeter.png'));
    imwrite(imoverlay(im_BW,im_P,'red'),strcat(OverlayPath,ImageName,'_BW_OverlayPerimeter.png'));
    
    % GT
    im_GT = imoverlay_v2(im,GT,'color',[0 1 0],'alpha',0.75,'bright');
    imwrite( im_GT, strcat(OverlayPath,ImageName,'_AN_GT_Overlay.png') );
    im_P = imdilate(boundarymask(GT), ones(3,3));
    imwrite(imoverlay(im,im_P,'green'),strcat(OverlayPath,ImageName,'_AN_GT_Perimeter.png'));
    imwrite(imoverlay(im_GT,im_P,'red'),strcat(OverlayPath,ImageName,'_AN_GT_OverlayPerimeter.png'));
    
    [OL_RGB,OL_IM] = FuseGTAndBW(GT,BW,im);
    imwrite(OL_RGB,strcat(OverlayPath,ImageName,'_AN_Overlay_BWonGT_RGB.png'));
    imwrite(OL_IM,strcat(OverlayPath,ImageName,'_AN_Overlay_BWonGT.png'));

    OL_RGB_OL = imoverlay_v2(im,OL_RGB(:,:,1),'color',[0 1 1],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,2),'color',[0 1 0],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,3),'color',[1 1 0],'alpha',0.75,'bright');
    imwrite(OL_RGB_OL,strcat(OverlayPath,ImageName,'_AN_Overlay_BWonGT_TPFPFN.png'));
    
%% Case 2: Strong Nuclei
    GT_Nuclei = GetOneClass(GT, ClassLabels, 'nucleus');
    GT_F = regionprops(GT_Nuclei,'Area','Centroid','BoundingBox',...
                                            'PixelIdxList','PixelList');
    N = length(GT_F);
    NR = CreateObjectEvaluationTable(N);
    
    % A array of (intially) zeros for correcting detecting segmented nuclei 
    BW_NucleiFound = zeros(length(BW_F),1);
    
    % loop for each nuclei in ground truth mask
    for j = 1:N
        found = false;
        % loop for each nuclei in segmentet mask
        for k=1:length(BW_F)
            % if the centroids of segmented nuclei inside GT nuclei
            if(sum(GT_F(j).PixelList(:,1) == int16(BW_F(k).Centroid(1))) > 0 && ...
                   sum(GT_F(j).PixelList(:,2) == int16(BW_F(k).Centroid(2))) > 0 && ...
                   BW_NucleiFound(k) == 0)
                BW_NucleiFound(k) = 1;
                found = true;
                NR.N_Found(j) = 1;
                IR_SN.D_TP(i) = IR_SN.D_TP(i) + 1;

                % Segmentation Metrics at Nuclei Level
                [NR.S_TP(j),NR.S_FN(j),NR.S_FP(j),NR.S_TPR(j),NR.S_PPV(j), ...
                    NR.S_FM(j),NR.S_OL(j),NR.S_SC(j)] ... 
                    = SegEva_Object( GT_F(j).PixelIdxList, BW_F(k).PixelIdxList);

                NR.S_Dic(j) = 2 * NR.S_TP(j) / ( GT_F(j).Area + BW_F(k).Area );

                % Segmentation Metric at Image Level
                alpha = GT_F(j).Area / sum(cat(1,GT_F.Area));
                beta = BW_F(k).Area / sum(cat(1,BW_F.Area));
                IR_SN.S_Dic_N(i) = IR_SN.S_Dic_N(i) + ...
                        ( ( (alpha * NR.S_Dic(j)) + (beta * NR.S_Dic(j)) ) / 2 );

                break;
            end
        end
        if(~found)
            IR_SN.D_FN(i) = IR_SN.D_FN(i) + 1;
            NR.N_Found(j) = -1;
        end
    end
    
    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( ResultPath, ImageName, '_Metrics_StrongNuclei.csv'));
    
    % Detection Metrics at Image Level
    IR_SN.D_FP(i) = length(BW_F) - IR_SN.D_TP(i);
    IR_SN.D_TPR(i) = IR_SN.D_TP(i) / (IR_SN.D_TP(i) + IR_SN.D_FN(i));
    IR_SN.D_PPV(i) = IR_SN.D_TP(i) / (IR_SN.D_TP(i) + IR_SN.D_FP(i));
    IR_SN.D_FM(i) = 2 * IR_SN.D_TPR(i) * IR_SN.D_PPV(i) / ...
                        (IR_SN.D_TPR(i) + IR_SN.D_PPV(i));

    % Segmentation Metrics at Image Level
    [ IR_SN.S_TPR(i),IR_SN.S_PPV(i),IR_SN.S_FM(i),IR_SN.S_TNR(i),IR_SN.S_FPR(i), ...
        IR_SN.S_AUC(i),IR_SN.S_Acc(i),IR_SN.S_OL(i),IR_SN.S_Dic(i),IR_SN.S_Kappa(i), ...
        IR_SN.S_SC(i),IR_SN.S_RI(i),IR_SN.S_ARI(i),IR_SN.S_GCE(i), ...
        IR_SN.S_MI(i),IR_SN.S_VI(i) ] = SegEva_Image(GT_Nuclei,BW);


    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( ResultPath, ImageName, '_Metrics_StrongNuclei.csv'));
    imwrite(GT_Nuclei, strcat(OverlayPath,ImageName,'_StrongNuclei.png'));

    % GT_StrongNuclei
    im_GT = imoverlay_v2(im,GT_Nuclei,'color',[0 1 0],'alpha',0.75,'bright');
    imwrite( im_GT, strcat(OverlayPath,ImageName,'_SN_GT_Overlay.png') );
    im_P = imdilate(boundarymask(GT_Nuclei), ones(3,3));
    imwrite(imoverlay(im,im_P,'green'),strcat(OverlayPath,ImageName,'_SN_GT_Perimeter.png'));
    imwrite(imoverlay(im_GT,im_P,'red'),strcat(OverlayPath,ImageName,'_SN_GT_OverlayPerimeter.png'));

    [OL_RGB,OL_IM] = FuseGTAndBW(GT_Nuclei,BW,im);
    imwrite(OL_RGB,strcat(OverlayPath,ImageName,'_SN_Overlay_BWonGT_RGB.png'));
    imwrite(OL_IM,strcat(OverlayPath,ImageName,'_SN_Overlay_BWonGT.png'));

    OL_RGB_OL = imoverlay_v2(im,OL_RGB(:,:,1),'color',[0 1 1],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,2),'color',[0 1 0],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,3),'color',[1 1 0],'alpha',0.75,'bright');
    imwrite(OL_RGB_OL,strcat(OverlayPath,ImageName,'_SN_Overlay_BWonGT_TPFPFN.png'));
    
%% Case 3; Weak Nuclei
    GT_Nuclei = GetOneClass(GT, ClassLabels, 'dimnucleus');
    GT_F = regionprops(GT_Nuclei,'Area','Centroid','BoundingBox',...
                                            'PixelIdxList','PixelList');
    N = length(GT_F);
    NR = CreateObjectEvaluationTable(N);
    
    % A array of (intially) zeros for correcting detecting segmented nuclei 
    BW_NucleiFound = zeros(length(BW_F),1);
    
    % loop for each nuclei in ground truth mask
    for j = 1:N
        found = false;
        % loop for each nuclei in segmentet mask
        for k=1:length(BW_F)
            % if the centroids of segmented nuclei inside GT nuclei
            if(sum(GT_F(j).PixelList(:,1) == int16(BW_F(k).Centroid(1))) > 0 && ...
                   sum(GT_F(j).PixelList(:,2) == int16(BW_F(k).Centroid(2))) > 0 && ...
                   BW_NucleiFound(k) == 0)
                BW_NucleiFound(k) = 1;
                found = true;
                NR.N_Found(j) = 1;
                IR_WN.D_TP(i) = IR_WN.D_TP(i) + 1;

                % Segmentation Metrics at Nuclei Level
                [NR.S_TP(j),NR.S_FN(j),NR.S_FP(j),NR.S_TPR(j),NR.S_PPV(j), ...
                    NR.S_FM(j),NR.S_OL(j),NR.S_SC(j)] ... 
                    = SegEva_Object( GT_F(j).PixelIdxList, BW_F(k).PixelIdxList);

                NR.S_Dic(j) = 2 * NR.S_TP(j) / ( GT_F(j).Area + BW_F(k).Area );

                % Segmentation Metric at Image Level
                alpha = GT_F(j).Area / sum(cat(1,GT_F.Area));
                beta = BW_F(k).Area / sum(cat(1,BW_F.Area));
                IR_WN.S_Dic_N(i) = IR_WN.S_Dic_N(i) + ...
                        ( ( (alpha * NR.S_Dic(j)) + (beta * NR.S_Dic(j)) ) / 2 );

                break;
            end
        end
        if(~found)
            IR_WN.D_FN(i) = IR_WN.D_FN(i) + 1;
            NR.N_Found(j) = -1;
        end
    end
    
    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( ResultPath, ImageName, '_Metrics_WeakNuclei.csv'));
    
    % Detection Metrics at Image Level
    IR_WN.D_FP(i) = length(BW_F) - IR_WN.D_TP(i);
    IR_WN.D_TPR(i) = IR_WN.D_TP(i) / (IR_WN.D_TP(i) + IR_WN.D_FN(i));
    IR_WN.D_PPV(i) = IR_WN.D_TP(i) / (IR_WN.D_TP(i) + IR_WN.D_FP(i));
    IR_WN.D_FM(i) = 2 * IR_WN.D_TPR(i) * IR_WN.D_PPV(i) / ...
                        (IR_WN.D_TPR(i) + IR_WN.D_PPV(i));

    % Segmentation Metrics at Image Level
    [ IR_WN.S_TPR(i),IR_WN.S_PPV(i),IR_WN.S_FM(i),IR_WN.S_TNR(i),IR_WN.S_FPR(i), ...
        IR_WN.S_AUC(i),IR_WN.S_Acc(i),IR_WN.S_OL(i),IR_WN.S_Dic(i),IR_WN.S_Kappa(i), ...
        IR_WN.S_SC(i),IR_WN.S_RI(i),IR_WN.S_ARI(i),IR_WN.S_GCE(i), ...
        IR_WN.S_MI(i),IR_WN.S_VI(i) ] = SegEva_Image(GT_Nuclei,BW);

    % Saving each Nuclei results (Nuclei level Results)
    writetable(NR,strcat( ResultPath, ImageName, '_Metrics_WeakNuclei.csv'));
    imwrite(GT_Nuclei, strcat(OverlayPath,ImageName,'_WeakNuclei.png'));
    
    % GT_WeakNuclei
    im_GT = imoverlay_v2(im,GT_Nuclei,'color',[0 1 0],'alpha',0.75,'bright');
    imwrite( im_GT, strcat(OverlayPath,ImageName,'_WN_GT_Overlay.png') );
    im_P = imdilate(boundarymask(GT_Nuclei), ones(3,3));
    imwrite(imoverlay(im,im_P,'green'),strcat(OverlayPath,ImageName,'_WN_GT_Perimeter.png'));
    imwrite(imoverlay(im_GT,im_P,'red'),strcat(OverlayPath,ImageName,'_WN_GT_OverlayPerimeter.png'));

    [OL_RGB,OL_IM] = FuseGTAndBW(GT_Nuclei,BW,im);
    imwrite(OL_RGB,strcat(OverlayPath,ImageName,'_WN_Overlay_BWonGT_RGB.png'));
    imwrite(OL_IM,strcat(OverlayPath,ImageName,'_WN_Overlay_BWonGT.png'));

    OL_RGB_OL = imoverlay_v2(im,OL_RGB(:,:,1),'color',[0 1 1],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,2),'color',[0 1 0],'alpha',0.75,'bright');
    OL_RGB_OL = imoverlay_v2(OL_RGB_OL,OL_RGB(:,:,3),'color',[1 1 0],'alpha',0.75,'bright');
    imwrite(OL_RGB_OL,strcat(OverlayPath,ImageName,'_WN_Overlay_BWonGT_TPFPFN.png'));    
end

% Saving each image results (Image Level Metrics)
writetable(IR,strcat( ResultPath, 'Image_DetectionAndSegmentationMetrics_AllNuclei.csv'));
writetable(IR_SN,strcat( ResultPath, 'Image_DetectionAndSegmentationMetrics_StrongNuclei.csv'));
writetable(IR_WN,strcat( ResultPath, 'Image_DetectionAndSegmentationMetrics_WeakNuclei.csv'));

    