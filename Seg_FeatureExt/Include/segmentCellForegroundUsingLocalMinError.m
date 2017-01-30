function [imCellForegroundMask] = segmentCellForegroundUsingLocalMinError(imInput, localWindowRadius, varargin)

    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])));
    p.addRequired('localWindowRadius', @(x) isscalar(x) && isnumeric(x));
    p.parse(imInput, localWindowRadius);

    p.addParameter('model', 'poisson', @(x) (ismember(x, {'gaussian', 'poisson'})));
    p.addParameter('localWindowPace', round(localWindowRadius / 3), @(x) isscalar(x) && isnumeric(x));
    p.addParameter('minLocalGlobalThresholdRatio', 0.5, @(x) isscalar(x) && isnumeric(x));    
    p.addParameter('minSliceToStackThresholdRatio', 0.4, @(x) isscalar(x) && isnumeric(x));
    p.addParameter('numHistogramBins', 256, @(x) (isscalar(x) && (x-floor(x)) == 0));
    p.addParameter('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
    p.addParameter('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
    p.addParameter('poolSize', 4, @(x) isscalar(x) && isnumeric(x)); 
    p.parse(imInput, localWindowRadius, varargin{:});
    
    strPdfModel = p.Results.model;    
    localWindowPace = round(p.Results.localWindowPace);
    minSliceLocalGlobalThresholdRatio = p.Results.minLocalGlobalThresholdRatio;   
    numHistogramBins = p.Results.numHistogramBins;
    minSliceToStackThresholdRatio = p.Results.minSliceToStackThresholdRatio;
    flagDebugMode = p.Results.debugMode;
    flagParallelize = p.Results.flagParallelize;
    poolSize = p.Results.poolSize;

    % local otsu thresholding in each slice
    imCellForegroundMask = zeros(size(imInput));
    
    threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                'numHistogramBins', numHistogramBins));
    globalStackThresh = threshFunc( imInput );
    
    if flagDebugMode
        imLocalThresholdVals = zeros(size(imInput));
        imGlobalSliceThresholdVals = zeros(size(imInput));
    end
    
    % perform thresholding
    if flagParallelize

        % if pool not open, then open pool
        if isempty(gcp('nocreate')) 
            poolobj = parpool(poolSize);
        end
        
        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Local Min-error Thresholding with a %s model in each of the %d slices ... ', strPdfModel, size(imInput, 3));
        end                                            
        
        parfor sliceId = 1:size(imInput, 3)
            
            imSlice = imInput(:,:,sliceId);

            threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                        'numHistogramBins', numHistogramBins));

            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minSliceToStackThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
                fprintf('\nWARNING: Slice - %d: SliceThreshold/StackThreshold of %.2f is less than the allowed lower-bound of %.2f ... \n', ...
                         sliceId, globalThreshVal/globalStackThresh, minSliceToStackThresholdRatio);
            end
                                                        
            imCellForegroundMask(:,:,sliceId) = imMask;
            
        end

        if flagDebugMode
            timeElapsed = toc(totalTimer);
            fprintf('took %f seconds\n', timeElapsed);
        end            
        
        % close matlab pool
         delete(poolobj);
        
    else
        
        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Local Min-error Thresholding with a %s model in each of the %d slices ... \n', strPdfModel, size(imInput, 3));
        end
        
        for sliceId = 1:size(imInput, 3)

            if flagDebugMode
                tic
                fprintf('\n\tThresholding slice %d/%d ... ', sliceId, size(imInput, 3));
            end
            
            imSlice = imInput(:,:,sliceId);

            threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                        'numHistogramBins', numHistogramBins));

            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minSliceToStackThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
                fprintf('\nWARNING: Slice - %d: SliceThreshold/StackThreshold of %.2f is less than the allowed lower-bound of %.2f ... \n', ...
                         sliceId, globalThreshVal/globalStackThresh, minSliceToStackThresholdRatio);
            end
                                                        
            imCellForegroundMask(:,:,sliceId) = imMask;
            
            if flagDebugMode
                timeElapsed = toc;
                fprintf('took %f seconds\n', timeElapsed);
            end            
            
        end

        if flagDebugMode
            timeElapsed = toc(totalTimer);
            fprintf('\nTotal time taken - %f seconds\n', timeElapsed);
        end            
        
    end
    
    % post-processing
    diskRad = ones(1,ndims(imInput));
    diskRad(1:2) = 3;
    imCellForegroundMask = imopen(imCellForegroundMask, streldisknd(diskRad));

        
end