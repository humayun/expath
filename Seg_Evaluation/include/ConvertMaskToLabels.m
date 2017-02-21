function [ Label ] = ConvertMaskToLabels( mask )
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    [X,Y,N] = size(mask);
    Label = zeros(X,Y);
    for j = 1:N
        Nuclei = logical(mask(:,:,j));
        Label(Nuclei) = j;
    end
end

