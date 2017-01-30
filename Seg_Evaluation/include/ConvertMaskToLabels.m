function [ Label ] = ConvertMaskToLabels( mask )

    [X,Y,N] = size(mask);
    Label = zeros(X,Y);
    for j = 1:N
        Nuclei = logical(mask(:,:,j));
        Label(Nuclei) = j;
    end
end

