function [ RLFeatures ] = ComputeRLFeatures( RGB, BW, NoOfChannels, GrayLevels, FeaturesPath, ImageName)
%ComputeRLFeatures compute Run-Length features in seven color 
%channels (Red, Green, Blue, HSV(V), Lab(L), H&E(H) and BR
%
    if(nargin < 4)
        GrayLevels = 64;
    end
    
    [x,y] = size(BW);
    [L,num] = bwlabel(BW);
    bb = regionprops(L,'BoundingBox'); 
    bb = floor(cat(1,bb.BoundingBox));

    RGB = uint8(256*mat2gray(RGB));
    GS_or_RGB = size(RGB,3);
    
    if (GS_or_RGB == 1)
        RLFeatures = struct( ...
        'SRE_GS',0, 'LRE_GS',0, 'GLN_GS',0, 'RLN_GS',0, 'RP_GS',0, 'LGRE_GS',0,...
        'HGRE_GS',0, 'SRLGE_GS',0, 'SRHGE_GS',0, 'LRLGE_GS',0, 'LRHGE_GS',0);

        for i = 1:num
            obj=RGB;
            obj(~(L==i)) = NaN;
            top_x = bb(i,2);
            if(top_x<1)
                top_x = 1;
            end
            top_y = bb(i,1);
            if(top_y<1)
                top_y = 1;
            end
            width = top_x+bb(i,4);
            if(width > x)
                width = x;
            end
            height = top_y+bb(i,3);
            if(height > y)
                height = y;
            end
            patch = obj(top_x:width,top_y:height);
            [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                GrayLevels,'GrayLimits',[]); 
            stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
            RLFeatures.SRE_GS(i,1) = mean(stats(:,1));
            RLFeatures.LRE_GS(i,1) = mean(stats(:,2));
            RLFeatures.GLN_GS(i,1) = mean(stats(:,3));
            RLFeatures.RLN_GS(i,1) = mean(stats(:,4));
            RLFeatures.RP_GS(i,1) = mean(stats(:,5));
            RLFeatures.LGRE_GS(i,1) = mean(stats(:,6));
            RLFeatures.HGRE_GS(i,1)= mean(stats(:,7));
            RLFeatures.SRLGE_GS(i,1) = mean(stats(:,8));
            RLFeatures.SRHGE_GS(i,1) = mean(stats(:,9));
            RLFeatures.LRLGE_GS(i,1) = mean(stats(:,10));
            RLFeatures.LRHGE_GS(i,1) = mean(stats(:,11));
        end
        
     elseif(GS_or_RGB == 3)
        if(NoOfChannels == 4 && isnumeric(NoOfChannels))
         
             RLFeatures = struct( ...
              'SRE_R',0, 'LRE_R',0, 'GLN_R',0, 'RLN_R',0, 'RP_R',0, 'LGRE_R',0, ...
              'HGRE_R',0, 'SRLGE_R',0, 'SRHGE_R',0, 'LRLGE_R',0, 'LRHGE_R',0, ...
              'SRE_V',0, 'LRE_V',0, 'GLN_V',0, 'RLN_V',0, 'RP_V',0, 'LGRE_V',0, ...
              'HGRE_V',0, 'SRLGE_V',0, 'SRHGE_V',0, 'LRLGE_V',0, 'LRHGE_V',0, ...
              'SRE_L',0, 'LRE_L',0, 'GLN_L',0, 'RLN_L',0, 'RP_L',0, 'LGRE_L',0, ...
              'HGRE_L',0, 'SRLGE_L',0, 'SRHGE_L',0, 'LRLGE_L',0, 'LRHGE_L',0, ...
              'SRE_H',0, 'LRE_H',0, 'GLN_H',0, 'RLN_H',0, 'RP_H',0, 'LGRE_H',0, ...
              'HGRE_H',0, 'SRLGE_H',0, 'SRHGE_H',0, 'LRLGE_H',0, 'LRHGE_H',0);

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_R(i,1) = mean(stats(:,1));
                RLFeatures.LRE_R(i,1) = mean(stats(:,2));
                RLFeatures.GLN_R(i,1) = mean(stats(:,3));
                RLFeatures.RLN_R(i,1) = mean(stats(:,4));
                RLFeatures.RP_R(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_R(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_R(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_R(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_R(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_R(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_R(i,1) = mean(stats(:,11));
            end

            %% 2 - V (HSV) Channel
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_V(i,1) = mean(stats(:,1));
                RLFeatures.LRE_V(i,1) = mean(stats(:,2));
                RLFeatures.GLN_V(i,1) = mean(stats(:,3));
                RLFeatures.RLN_V(i,1) = mean(stats(:,4));
                RLFeatures.RP_V(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_V(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_V(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_V(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_V(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_V(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_V(i,1) = mean(stats(:,11));
            end

            %% 3 - L (Lab) Channel
            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_L(i,1) = mean(stats(:,1));
                RLFeatures.LRE_L(i,1) = mean(stats(:,2));
                RLFeatures.GLN_L(i,1) = mean(stats(:,3));
                RLFeatures.RLN_L(i,1) = mean(stats(:,4));
                RLFeatures.RP_L(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_L(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_L(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_L(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_L(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_L(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_L(i,1) = mean(stats(:,11));
            end

            %% 4 - H (H&E) Channel
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_H(i,1) = mean(stats(:,1));
                RLFeatures.LRE_H(i,1) = mean(stats(:,2));
                RLFeatures.GLN_H(i,1) = mean(stats(:,3));
                RLFeatures.RLN_H(i,1) = mean(stats(:,4));
                RLFeatures.RP_H(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_H(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_H(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_H(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_H(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_H(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_H(i,1) = mean(stats(:,11));
            end

        elseif(NoOfChannels == 7 && isnumeric(NoOfChannels))
             RLFeatures = struct( ...
              'SRE_R',0, 'LRE_R',0, 'GLN_R',0, 'RLN_R',0, 'RP_R',0, 'LGRE_R',0, ...
              'HGRE_R',0, 'SRLGE_R',0, 'SRHGE_R',0, 'LRLGE_R',0, 'LRHGE_R',0, ...
              'SRE_G',0, 'LRE_G',0, 'GLN_G',0, 'RLN_G',0, 'RP_G',0, 'LGRE_G',0, ...
              'HGRE_G',0, 'SRLGE_G',0, 'SRHGE_G',0, 'LRLGE_G',0, 'LRHGE_G',0, ...
              'SRE_B',0, 'LRE_B',0, 'GLN_B',0, 'RLN_B',0, 'RP_B',0, 'LGRE_B',0, ...
              'HGRE_B',0, 'SRLGE_B',0, 'SRHGE_B',0, 'LRLGE_B',0, 'LRHGE_B',0, ...
              'SRE_V',0, 'LRE_V',0, 'GLN_V',0, 'RLN_V',0, 'RP_V',0, 'LGRE_V',0, ...
              'HGRE_V',0, 'SRLGE_V',0, 'SRHGE_V',0, 'LRLGE_V',0, 'LRHGE_V',0, ...
              'SRE_L',0, 'LRE_L',0, 'GLN_L',0, 'RLN_L',0, 'RP_L',0, 'LGRE_L',0, ...
              'HGRE_L',0, 'SRLGE_L',0, 'SRHGE_L',0, 'LRLGE_L',0, 'LRHGE_L',0, ...
              'SRE_H',0, 'LRE_H',0, 'GLN_H',0, 'RLN_H',0, 'RP_H',0, 'LGRE_H',0, ...
              'HGRE_H',0, 'SRLGE_H',0, 'SRHGE_H',0, 'LRLGE_H',0, 'LRHGE_H',0, ...
              'SRE_Br',0,'LRE_Br',0,'GLN_Br',0,'RLN_Br',0,'RP_Br',0,'LGRE_Br',0,...
              'HGRE_Br',0, 'SRLGE_Br',0, 'SRHGE_Br',0, 'LRLGE_Br',0, 'LRHGE_Br',0);

            %% 1 - Red Channel
            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_R(i,1) = mean(stats(:,1));
                RLFeatures.LRE_R(i,1) = mean(stats(:,2));
                RLFeatures.GLN_R(i,1) = mean(stats(:,3));
                RLFeatures.RLN_R(i,1) = mean(stats(:,4));
                RLFeatures.RP_R(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_R(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_R(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_R(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_R(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_R(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_R(i,1) = mean(stats(:,11));
            end

            %% 2 - Green Channel
            channel = RGB(:,:,2);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_G(i,1) = mean(stats(:,1));
                RLFeatures.LRE_G(i,1) = mean(stats(:,2));
                RLFeatures.GLN_G(i,1) = mean(stats(:,3));
                RLFeatures.RLN_G(i,1) = mean(stats(:,4));
                RLFeatures.RP_G(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_G(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_G(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_G(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_G(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_G(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_G(i,1) = mean(stats(:,11));
            end

            %% 3 - Blue Channel
            channel = RGB(:,:,3);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_B(i,1) = mean(stats(:,1));
                RLFeatures.LRE_B(i,1) = mean(stats(:,2));
                RLFeatures.GLN_B(i,1) = mean(stats(:,3));
                RLFeatures.RLN_B(i,1) = mean(stats(:,4));
                RLFeatures.RP_B(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_B(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_B(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_B(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_B(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_B(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_B(i,1) = mean(stats(:,11));
            end

            %% 4 - V (HSV) Channel
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_V(i,1) = mean(stats(:,1));
                RLFeatures.LRE_V(i,1) = mean(stats(:,2));
                RLFeatures.GLN_V(i,1) = mean(stats(:,3));
                RLFeatures.RLN_V(i,1) = mean(stats(:,4));
                RLFeatures.RP_V(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_V(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_V(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_V(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_V(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_V(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_V(i,1) = mean(stats(:,11));
            end

            %% 5 - L (Lab) Channel
            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_L(i,1) = mean(stats(:,1));
                RLFeatures.LRE_L(i,1) = mean(stats(:,2));
                RLFeatures.GLN_L(i,1) = mean(stats(:,3));
                RLFeatures.RLN_L(i,1) = mean(stats(:,4));
                RLFeatures.RP_L(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_L(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_L(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_L(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_L(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_L(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_L(i,1) = mean(stats(:,11));
            end

            %% 6 - H (H&E) Channel
            %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
            %[ DCH, H, E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_H(i,1) = mean(stats(:,1));
                RLFeatures.LRE_H(i,1) = mean(stats(:,2));
                RLFeatures.GLN_H(i,1) = mean(stats(:,3));
                RLFeatures.RLN_H(i,1) = mean(stats(:,4));
                RLFeatures.RP_H(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_H(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_H(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_H(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_H(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_H(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_H(i,1) = mean(stats(:,11));
            end

            %% 7 - BlueRatio Image
            channel= uint8(GrayLevels*mat2gray(ComputeBlueRatio(RGB)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_Br(i,1) = mean(stats(:,1));
                RLFeatures.LRE_Br(i,1) = mean(stats(:,2));
                RLFeatures.GLN_Br(i,1) = mean(stats(:,3));
                RLFeatures.RLN_Br(i,1) = mean(stats(:,4));
                RLFeatures.RP_Br(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_Br(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_Br(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_Br(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_Br(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_Br(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_Br(i,1) = mean(stats(:,11));
            end
            
        elseif(NoOfChannels == 'R' && ~isnumeric(NoOfChannels))
            
            %% 1 - Red Channel
             RLFeatures = struct( ...
              'SRE_R',0, 'LRE_R',0, 'GLN_R',0, 'RLN_R',0, 'RP_R',0, 'LGRE_R',0, ...
              'HGRE_R',0, 'SRLGE_R',0, 'SRHGE_R',0, 'LRLGE_R',0, 'LRHGE_R',0);

            channel = RGB(:,:,1);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_R(i,1) = mean(stats(:,1));
                RLFeatures.LRE_R(i,1) = mean(stats(:,2));
                RLFeatures.GLN_R(i,1) = mean(stats(:,3));
                RLFeatures.RLN_R(i,1) = mean(stats(:,4));
                RLFeatures.RP_R(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_R(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_R(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_R(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_R(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_R(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_R(i,1) = mean(stats(:,11));
            end
        elseif(NoOfChannels == 'G' && ~isnumeric(NoOfChannels))
            
            %% Green Channel
             RLFeatures = struct( ...
              'SRE_G',0, 'LRE_G',0, 'GLN_G',0, 'RLN_G',0, 'RP_G',0, 'LGRE_G',0, ...
              'HGRE_G',0, 'SRLGE_G',0, 'SRHGE_G',0, 'LRLGE_G',0, 'LRHGE_G',0);

            channel = RGB(:,:,2);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_G(i,1) = mean(stats(:,1));
                RLFeatures.LRE_G(i,1) = mean(stats(:,2));
                RLFeatures.GLN_G(i,1) = mean(stats(:,3));
                RLFeatures.RLN_G(i,1) = mean(stats(:,4));
                RLFeatures.RP_G(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_G(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_G(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_G(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_G(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_G(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_G(i,1) = mean(stats(:,11));
            end
            
        elseif(NoOfChannels == 'B' && ~isnumeric(NoOfChannels))
            
            %% Blue Channel
             RLFeatures = struct( ...
              'SRE_B',0, 'LRE_B',0, 'GLN_B',0, 'RLN_B',0, 'RP_B',0, 'LGRE_B',0, ...
              'HGRE_B',0, 'SRLGE_B',0, 'SRHGE_B',0, 'LRLGE_B',0, 'LRHGE_B',0);

            [x,y] = size(BW);
            [L,num] = bwlabel(BW);
            bb = regionprops(L,'BoundingBox'); 
            bb = floor(cat(1,bb.BoundingBox));
            
            channel = RGB(:,:,3);
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_B(i,1) = mean(stats(:,1));
                RLFeatures.LRE_B(i,1) = mean(stats(:,2));
                RLFeatures.GLN_B(i,1) = mean(stats(:,3));
                RLFeatures.RLN_B(i,1) = mean(stats(:,4));
                RLFeatures.RP_B(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_B(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_B(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_B(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_B(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_B(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_B(i,1) = mean(stats(:,11));
            end

        elseif(NoOfChannels == 'V' && ~isnumeric(NoOfChannels))
            
            %% V (HSV) Channel
             RLFeatures = struct( ...
              'SRE_V',0, 'LRE_V',0, 'GLN_V',0, 'RLN_V',0, 'RP_V',0, 'LGRE_V',0, ...
              'HGRE_V',0, 'SRLGE_V',0, 'SRHGE_V',0, 'LRLGE_V',0, 'LRHGE_V',0);
          
            channel = rgb2hsv(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,3)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_V(i,1) = mean(stats(:,1));
                RLFeatures.LRE_V(i,1) = mean(stats(:,2));
                RLFeatures.GLN_V(i,1) = mean(stats(:,3));
                RLFeatures.RLN_V(i,1) = mean(stats(:,4));
                RLFeatures.RP_V(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_V(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_V(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_V(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_V(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_V(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_V(i,1) = mean(stats(:,11));
            end
            
        elseif(NoOfChannels == 'L' && ~isnumeric(NoOfChannels))
            
            %% L (Lab) Channel
             RLFeatures = struct( ...
              'SRE_L',0, 'LRE_L',0, 'GLN_L',0, 'RLN_L',0, 'RP_L',0, 'LGRE_L',0, ...
              'HGRE_L',0, 'SRLGE_L',0, 'SRHGE_L',0, 'LRLGE_L',0, 'LRHGE_L',0);
            
            channel = rgb2lab(RGB);
            channel = uint8(GrayLevels*mat2gray(channel(:,:,1)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_L(i,1) = mean(stats(:,1));
                RLFeatures.LRE_L(i,1) = mean(stats(:,2));
                RLFeatures.GLN_L(i,1) = mean(stats(:,3));
                RLFeatures.RLN_L(i,1) = mean(stats(:,4));
                RLFeatures.RP_L(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_L(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_L(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_L(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_L(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_L(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_L(i,1) = mean(stats(:,11));
            end
            
        elseif(NoOfChannels == 'H' && ~isnumeric(NoOfChannels))
            
            %% H (H&E) Channel
             RLFeatures = struct( ...
              'SRE_H',0, 'LRE_H',0, 'GLN_H',0, 'RLN_H',0, 'RP_H',0, 'LGRE_H',0, ...
              'HGRE_H',0, 'SRLGE_H',0, 'SRHGE_H',0, 'LRLGE_H',0, 'LRHGE_H',0);

            [ channel, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
            channel = uint8(GrayLevels*mat2gray(channel));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_H(i,1) = mean(stats(:,1));
                RLFeatures.LRE_H(i,1) = mean(stats(:,2));
                RLFeatures.GLN_H(i,1) = mean(stats(:,3));
                RLFeatures.RLN_H(i,1) = mean(stats(:,4));
                RLFeatures.RP_H(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_H(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_H(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_H(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_H(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_H(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_H(i,1) = mean(stats(:,11));
            end
            
        elseif(strcmp(NoOfChannels,'BR') && ~isnumeric(NoOfChannels))
            %% BR Channel
             RLFeatures = struct( ...
              'SRE_Br',0,'LRE_Br',0,'GLN_Br',0,'RLN_Br',0,'RP_Br',0,'LGRE_Br',0,...
              'HGRE_Br',0, 'SRLGE_Br',0, 'SRHGE_Br',0, 'LRLGE_Br',0, 'LRHGE_Br',0);

            channel= uint8(GrayLevels*mat2gray(ComputeBlueRatio(RGB)));
            for i = 1:num
                obj=channel;
                obj(~(L==i)) = NaN;
                top_x = bb(i,2);
                if(top_x<1)
                    top_x = 1;
                end
                top_y = bb(i,1);
                if(top_y<1)
                    top_y = 1;
                end
                width = top_x+bb(i,4);
                if(width > x)
                    width = x;
                end
                height = top_y+bb(i,3);
                if(height > y)
                    height = y;
                end
                patch = obj(top_x:width,top_y:height);
                [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
                    GrayLevels,'GrayLimits',[]); 
                stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
                    'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
                RLFeatures.SRE_Br(i,1) = mean(stats(:,1));
                RLFeatures.LRE_Br(i,1) = mean(stats(:,2));
                RLFeatures.GLN_Br(i,1) = mean(stats(:,3));
                RLFeatures.RLN_Br(i,1) = mean(stats(:,4));
                RLFeatures.RP_Br(i,1) = mean(stats(:,5));
                RLFeatures.LGRE_Br(i,1) = mean(stats(:,6));
                RLFeatures.HGRE_Br(i,1)= mean(stats(:,7));
                RLFeatures.SRLGE_Br(i,1) = mean(stats(:,8));
                RLFeatures.SRHGE_Br(i,1) = mean(stats(:,9));
                RLFeatures.LRLGE_Br(i,1) = mean(stats(:,10));
                RLFeatures.LRHGE_Br(i,1) = mean(stats(:,11));
            end
            
        end
    end
    if(nargin == 6)
         struct2csv(RLFeatures,strcat(FeaturesPath,ImageName,'_RLFeatures.csv'));
    end
    RLFeatures = struct2table(RLFeatures);
end

