function [ Out, im_rgb ] = FuseGTAndBW( GT, BW, im )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    GT = logical(GT);
    BW = logical(BW);
    [Y,X] = size(GT);
    Out = zeros(Y,X,3);
    C = ndims(im);
    if C == 2
        im_rgb = gray2rgb(im);
    else
        im_rgb = im;
    end
    for i = 1:Y
        for j = 1:X
            if GT(i,j)==1 && BW(i,j)==1   % TP 
                im_rgb(i,j,:) = [0,255,0];
                Out(i,j,:) = [0,255,0];
            elseif GT(i,j)==1 && BW(i,j)==0 % FN
                im_rgb(i,j,:) = [255,0,0];
                Out(i,j,:) = [255,0,0];
            elseif GT(i,j)==0 && BW(i,j)==1 % FP
                im_rgb(i,j,:) = [0,0,255];
                Out(i,j,:) = [0,0,255];
            end
        end
    end
end

