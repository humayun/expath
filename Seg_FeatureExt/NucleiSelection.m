function [BW, NoSmallObjRemoved, NoLargeObjRemoved] = NucleiSelection(BW,MinSize,MaxSize,Conn)
%SIZEFILTERING function remove small and large object from binary image
%
% ########################################################################

    [~,L1] = bwlabeln(BW);
    BW = bwareaopen(BW,MinSize,Conn);
    [~,L2] = bwlabeln(BW);
    NoSmallObjRemoved = L1-L2;
    msg = ['Removed ' int2str(NoSmallObjRemoved) ' regions due to small size.'];
    disp(msg);
    % remove the large cells
    [~,L1] = bwlabeln(BW);
    BW = BW - bwareaopen(BW,MaxSize,Conn);
    [~,L2] = bwlabeln(BW);
    NoLargeObjRemoved = L1-L2;
    msg = ['Removed ' int2str(NoLargeObjRemoved) ' regions due to large size.'];
    disp(msg);
end