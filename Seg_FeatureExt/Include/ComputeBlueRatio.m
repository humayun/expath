function [ BR ] = ComputeBlueRatio( RGB )
%ComputeBlueRatio compute BlueRatio of RGB image
%   This function computes blueratio image from RGB image
%
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    n1 = 100.*(double(RGB(:,:,3)));
    d1 = (1+(double(RGB(:,:,1))+double(RGB(:,:,2))));
    d2 = (1+(double(RGB(:,:,1))+double(RGB(:,:,2))+double(RGB(:,:,3))));
    BR= n1./d1 .* 256./d2 ;

end

