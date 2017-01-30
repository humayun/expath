function [ imNucleiSeedPoints ] = SeedsSelection( imNucleiSeedPoints, p )
%SeedsSelection function select seed point for nuclei segmentation
%   Detailed explanation goes here
     
    % Suppress all maxima with value less than 0.01% of the maximum value
    maxVal = max(imNucleiSeedPoints(:));
    imNucleiSeedPoints(imNucleiSeedPoints<=p.thresholdvalue*maxVal) = 0;

    % remove regions without a seed point or with weak seed responses
    L = bwlabeln( imNucleiSeedPoints );
    fgRegProps = regionprops(L, 'PixelIdxList' );
    flagPruneRegion = false( numel(fgRegProps), 1 );
    for j = 1:numel(fgRegProps)

        curRgnSeedPixInd = fgRegProps(j).PixelIdxList( imNucleiSeedPoints( fgRegProps(j).PixelIdxList ) > 0 );
        
        % zero-out regions without seed-point
        if isempty(curRgnSeedPixInd)
            flagPruneRegion(j) = true;
            continue;
        end
        
         % zero-out regions with all seed points with response below a speficied cutoff
         if max( p.imFilterResponse(curRgnSeedPixInd) ) < p.seedDetectionResponseCutoff
             flagPruneRegion(j) = true;             
             continue;
         end
 
         % zero-out regions with max distance map value below a specified cutoff
         if max(p.fgDistanceMap(fgRegProps(j).PixelIdxList)) < p.seedToBackgroundDistanceCutoff
             flagPruneRegion(j) = true;             
             continue;
         end
    end
    L( ismember(L, find(flagPruneRegion)) ) = 0;
    imNucleiSeedPoints = double( L > 0 );

end

