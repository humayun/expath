function [ Mask, Count ] = GetOneClass( Mask, ClassLabels, ClassLabel )
%GetOneClass Summary of this function goes here
%   GetOneClass function take Mask image (2D), with class labels with each
%   object in it and selected class label and output the label of selected
%   class objects.
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    Count = 0;
    for i = 1:max(max(Mask))
        if( ~ strcmp(cellstr(ClassLabels{1,i}),ClassLabel))
            Mask (Mask == i) = 0;
        else
            Count = Count + 1;
        end
    end

end

