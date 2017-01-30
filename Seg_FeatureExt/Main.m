function  Main( InputFolder )
%Main is starting point of program
%   It will take input the folder address having tif files - 3D Lightsheet
%   microscopic images

	addpath('/home/hi41/WS/Exp_Segmentation_FeatureExtraction/Include')
	
    %% Inilialize Parameters
    p = InitializeParameters();

    %% Read 3D images
    SegPath = strcat(InputFolder,'Segmentation/');
    mkdir(SegPath);
    FeaturePath = strcat(InputFolder,'Features/');
    mkdir(FeaturePath);
    srcFiles = dir(strcat(InputFolder,'*.tif'));
    for i = 1:length(srcFiles)
        [~, ImageName, ~] = fileparts(srcFiles(i).name);
        im = imread(strcat(InputFolder,srcFiles(i).name));

        p.seedDetectionResponseCutoff = (4.0/4096) * max( im(:) );
        
        % Image Enhancement
        im = medfilt2(im,[10,10]);

        %% Foreground Segmentation
        %fg = im2bw(im, graythresh(im)*0.85);
        %fg = double(im > thresholdOtsu(im)*0.7);
        fg = segmentCellForegroundUsingLocalMinError(im, p.localWindowRadius, ...
            'model', p.threshmodel, ...
            'localWindowPace', p.localWindowPace, ...
            'minLocalGlobalThresholdRatio', p.minLocalGlobalThresholdRatio, ...
            'numHistogramBins', p.numHistogramBins, ...
            'debugMode', p.debugMode, ...
            'flagParallelize',p.flagParallelize, ...
            'poolSize', p.poolSize);

        % fill holes and remove small regions by morphological operations
        fg = imclose(fg, streldisknd(2*ones(1,ndims(im))) );

        % remove regions with small and invalid/unusual sizes
        L = bwlabeln( fg );
        stats = regionprops( L, {'Area', 'BoundingBox'} );
        validBBox = false(1, numel(stats));
        imsize = size(im);
        if(~isempty(stats))
            for j = 1:numel( stats )
                bboxSideLength = stats(j).BoundingBox((ndims(im)+1):end);
                bBoxBigEnough = all( bboxSideLength >= p.minNucleiBBox );
                if ~bBoxBigEnough
                    continue;
                end
                bboxSideLength = bboxSideLength(1:2);
                IsGridRegion = any(bboxSideLength >= 0.5 * imsize([2,1])) ...
                    && min(bboxSideLength) / max(bboxSideLength) < 0.25;
                validBBox(j) = ~IsGridRegion;
            end
            regArea = [stats.Area] * prod(p.spacing);
            smallRegInd = find( ~validBBox | regArea < p.minNucleiArea );
            L( ismember(L, smallRegInd) ) = 0;
            fg = double(L > 0);

            % removes thin vessel like structures
            % 0.20 (more smaller regions) or 0.3 (less smaller regions)
            cleanDiskRad=max([round(0.20*min(p.nucleiDiameterRange)./p.spacing); [2,2]] );
            fg = imopen(fg, streldisknd(cleanDiskRad) ); 

            p.fgDistanceMap = bwdistsc( 1-fg, p.spacing );
            p.numLoGScales = min(p.maxNumLoGScales, p.nucleiDiameterRange(2)-p.nucleiDiameterRange(1));

            %% Seed Detection
            [imNucleiSeedPoints,p.imFilterResponse,~] = ...
                    detectBlobsUsingAdaptiveMultiscaleLoG(im, p.fgDistanceMap, ...
                    'blobDiameterRange', p.nucleiDiameterRange, ...       
                    'flagParallelize', p.flagParallelize, ...
                    'minBlobDistance', p.minDistanceBWSeeds, ...
                    'logResponseCutoff',p.logResponseCutoff, ...
                    'numLoGScales', p.numLoGScales, ...
                    'poolSize', p.poolSize);

            imNucleiSeedPoints  = SeedsSelection(  imNucleiSeedPoints, p );

            %% Marker-Controlled Watershed based nuclei segmentation and separation

            imFilterResponse = 1 - mat2gray( p.imFilterResponse );
            imFilterResponse(~fg) = 100; % some high value
            imMinimaImposedAtSeeds = imimposemin( imFilterResponse, imNucleiSeedPoints );
            L = watershed( imMinimaImposedAtSeeds );
            L( ~fg ) = 0; % Remove background using foreground mask

            minObjDiameterImsp = min(p.nucleiDiameterRange(1) ./ p.spacing(1:2));    
            cleanStrel = strel('disk', round(0.25 * minObjDiameterImsp));
            L = imopen(L > 0, cleanStrel);       
            L = bwlabeln( L );

            regLabelWithSeed = unique( L( imNucleiSeedPoints > 0 ) );
            L( ~ismember(L, regLabelWithSeed) ) = 0; 
            L = (L > 0);

            [L,~,~] = NucleiSelection(L,p.minNucleiArea,p.maxNucleiArea,8);

            imwrite(L, strcat(SegPath, ImageName, '_Binary.png'));

            %% Create green color contour of nuclei
            [x,y] = size(im);
            overlay = zeros(x,y,3);
            overlay(:,:,2) = bwperim(L);
            overlay = imfuse(im, overlay, 'blend');
            imwrite(overlay, strcat(SegPath, ImageName, '_Overlay.png'));
            %figure, imshow(overlay)

            ComputeMorphologicalFeatures( L, FeaturePath, ImageName );
            ComputeIntensityFeatures(im, L, 0, FeaturePath, ImageName);
            ComputeCMFeatures( im, L, 0,  p.GrayLevels, FeaturePath, ImageName);
            ComputeRLFeatures( im, L, 0,  p.GrayLevels, FeaturePath, ImageName);

        end
    end
	
	%% Compute features summaries 
	ComputeFeaturesSummary(FeaturePath);
end

