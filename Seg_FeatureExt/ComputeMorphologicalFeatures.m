function [MorphFeatures, NucleiCentroids] = ComputeMorphologicalFeatures( BW, FeaturePath, ImageName )
%ComputeMorphologicalFeatures compute morphological features

    morphFeatures = regionprops(logical(BW),'Centroid', 'Area', ...
        'Perimeter', 'MajorAxisLength', 'MinorAxisLength', ...
        'Eccentricity', 'ConvexArea', 'Orientation', ...
        'EquivDiameter', 'Solidity', 'Extent', 'Perimeter', 'PixelIdxList');
    
    Centroid = cat(1,morphFeatures.Centroid);
    CentroidX = uint16(Centroid(:,1));
    CentroidY = uint16(Centroid(:,2));
    NucleiCentroids.Centroid_X = CentroidX;
    NucleiCentroids.Centroid_Y = CentroidY;

    MorphFeatures = struct('Area',0, 'Perimeter',0,'MajorAxisLength',0, ...
        'MinorAxisLength',0, 'Eccentricity',0, 'ConvexArea',0, ...
        'Orientation',0, 'EquivDiameter',0, 'Solidity',0, 'Extent',0, ...
        'Compactness',0, 'Ellipse_X',0, 'Ellipse_Y',0 );

    Area = cat(1,morphFeatures.Area);
    Perimeter = cat(1,morphFeatures.Perimeter);
    Compactness = 4*pi*Area./(Perimeter).^2;

    MorphFeatures.Area = Area;
    MorphFeatures.Perimeter = Perimeter;
    MorphFeatures.MajorAxisLength = cat(1,morphFeatures.MajorAxisLength);
    MorphFeatures.MinorAxisLength = cat(1,morphFeatures.MinorAxisLength);
    MorphFeatures.Eccentricity = cat(1,morphFeatures.Eccentricity);
    MorphFeatures.ConvexArea = cat(1,morphFeatures.ConvexArea);
    MorphFeatures.Orientation = cat(1,morphFeatures.Orientation);
    MorphFeatures.EquivDiameter = cat(1,morphFeatures.EquivDiameter);
    MorphFeatures.Solidity = cat(1,morphFeatures.Solidity);
    MorphFeatures.Extent = cat(1,morphFeatures.Extent);
    MorphFeatures.Compactness = Compactness;

    arrEllipse_X = zeros(numel(morphFeatures),1);
    arrEllipse_Y = zeros(numel(morphFeatures),1);
    
    % fit ellipse and store its radii
    for i = 1:numel(morphFeatures)

        [yind, xind] = ind2sub( size(BW), morphFeatures(i).PixelIdxList ); 
        ptCell = [xind, yind] .* repmat( size(BW), [numel(xind), 1] );
        ptCell = ptCell - repmat( mean(ptCell), [size(ptCell,1), 1] );
        [~, S, ~] = svd( (ptCell' * ptCell) / size(ptCell,1) );

         arrEllipse_X(i) = 2 * sqrt(S(1,1)); % eigen-values are a measure of variance
         arrEllipse_Y(i) = 2 * sqrt(S(2,2)); % eigen-values are a measure of variance
    end
    MorphFeatures.Ellipse_X = arrEllipse_X;
    MorphFeatures.Ellipse_Y = arrEllipse_Y;
    if(nargin == 3)
        struct2csv(NucleiCentroids,strcat(FeaturePath,ImageName,'_Centroid.csv'));
        struct2csv(MorphFeatures,strcat(FeaturePath,ImageName,'_MorphologicalFeatures.csv'));
    end
    MorphFeatures = struct2table(MorphFeatures);
end