function I = imoverlay_v2(varargin)
%IMOVERLAY Adds an overlay mask to a color image.
%
%   I1 = imoverlay(I, MASK, ...) adds an overlay mask MASK to an image I 
%   and outputs the resuting image I1.  MASK is typically a binary image, 
%   but an RGB color mask may also be used if desired.  Optional argument 
%   'color' allows the user to specify the color of the overlay mask 
%   ([0.4 0] by default).  The color must be specified as a RGB vector with
%   values between 0 and 1.  Optional argument 'alpha' (default 0.75) 
%   allows the user to specify the blending strength of the overlay, which 
%   ranges from weak (0) to strong (1). Optional argument 'interp' allows 
%   the user to specify the interpolation method used if MASK needs to be 
%   resized ('bicubic', 'nearest', etc).  The option 'bright' will enhance
%   the overlay mask for bright images.
%   
%   Example:
%   ------------
%   I = imread('cameraman.tif');
%   MASK = imcomplement(imread('testpat1.png'));
%   I1 = imoverlay(I, MASK, 'color', [ 0 .5 0], 'alpha', .8);
%   figure; imshow(I1);
%
%   Example 2:
%   ------------
%   I = imread('saturn.png');
%   MASK = imread('rice.png') > 150;
%   MASK = label2rgb(bwlabel(MASK), 'jet', [0 0 0]);
%   I1 = imoverlay(I, MASK, 'alpha', .5, 'interp', 'bicubic');
%   figure; imshow(I1);
%   
%   Copyright 2008 Kevin Smith
%
%   See also GRAY2RGB, IMOVERLAY, IMADD, IMLINCOMB, MIP, MAT2GRAY


I = varargin{1};
MASK = varargin{2};
color = [ 0 .4 0]; 
alpha = .75;
b_enhance = 0;

if nargin > 2;
    for i = 3:nargin
        if strcmp(varargin{i}, 'color')
            color = varargin{i+1};
        end
        if strcmp(varargin{i}, 'alpha')
            alpha = varargin{i+1};
        end
        if strcmp(varargin{i}, 'interp')
            interp = varargin{i+1};
        end
        if strcmp(varargin{i}, 'bright')
            b_enhance = 1;
        end
    end
end
        
% if I is not an RGB image, make it RGB
if size(I,3) ~= 3;
    I = gray2rgb(I);
end


% resize the mask to fit I, if necessary
if ~isequal(size(I,1), size(MASK,1))  || ~isequal(size(I,2), size(MASK,2)) 
    if exist('interp', 'var')
        MASK = imresize(MASK, [size(I,1) size(I,2)], interp);
    else
         MASK = imresize(MASK, [size(I,1) size(I,2)]);
    end
end


cls = class(I);                         % determine the data type of I


if length(size(MASK)) == 2
    % 2D MASK
    
    MASK2 = zeros([size(I,1) size(I,2) 3], cls);
    switch cls
        case 'uint8'
            MASK2(:,:,1) = uint8(color(1)*double(intmax(cls))*MASK);
            MASK2(:,:,2) = uint8(color(2)*double(intmax(cls))*MASK);
            MASK2(:,:,3) = uint8(color(3)*double(intmax(cls))*MASK);
            if b_enhance
                I(repmat(MASK, [1 1 3]) > 0) = uint8( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'uint16'
            MASK2(:,:,1) = uint16(color(1)*double(intmax(cls))*MASK);
            MASK2(:,:,2) = uint16(color(2)*double(intmax(cls))*MASK);
            MASK2(:,:,3) = uint16(color(3)*double(intmax(cls))*MASK);
            if b_enhance
                I(repmat(MASK, [1 1 3]) > 0) = uint16( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'uint32'
            MASK2(:,:,1) = uint32(color(1)*double(intmax(cls))*MASK);
            MASK2(:,:,2) = uint32(color(2)*double(intmax(cls))*MASK);
            MASK2(:,:,3) = uint32(color(3)*double(intmax(cls))*MASK);
            if b_enhance
                I(repmat(MASK, [1 1 3]) > 0) = uint32( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'double'
            MASK2(:,:,1) = color(1)*MASK;
            MASK2(:,:,2) = color(2)*MASK;
            MASK2(:,:,3) = color(3)*MASK;
            if b_enhance
                I(repmat(MASK, [1 1 3]) > 0) = (alpha)* I(repmat(MASK, [1 1 3]) > 0);
            end
    end
    
else
    % 3D MASK
    switch cls
        case 'uint8'
            MASK2 = uint8(MASK);
            if b_enhance
                MASK = max(MASK, [ ], 3);
                I(repmat(MASK, [1 1 3]) > 0) = uint8( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'uint16'
            MASK2 = uint16(MASK);
            if b_enhance
                MASK = max(MASK, [ ], 3);
                I(repmat(MASK, [1 1 3]) > 0) = uint16( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'uint32'
            MASK2 = uint32(MASK);
            if b_enhance
                MASK = max(MASK, [ ], 3);
                I(repmat(MASK, [1 1 3]) > 0) = uint32( (alpha)* I(repmat(MASK, [1 1 3]) > 0));
            end
        case 'double'
            MASK2 = double(MASK);
            if b_enhance
                MASK = max(MASK, [ ], 3);
                I(repmat(MASK, [1 1 3]) > 0) = (alpha)* I(repmat(MASK, [1 1 3]) > 0);
            end 
    end
end

I = imlincomb(1, I, alpha, MASK2, cls);