function GenerateMaskFromAnnotation()
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    %% Downloadiing User Database (Images + Annotations) from LabelMe Website
    % Define the root folder for the images
    HOMEIMAGES = 'PreExpanded'; 
    HOMEANNOTATIONS = 'PreExpanded/Annotation'; 
    HOMEMASKS = 'PreExpanded/Masks/'; 
    mkdir(HOMEMASKS);

    % Database name with path on LabelMe website
    DatabaseName ={'becklab_perexpanded'};

    LMinstall(DatabaseName, HOMEIMAGES,HOMEANNOTATIONS);
    
    %% Reading the index
    D = LMdatabase(HOMEANNOTATIONS);

    %% Extracting Mask for all objects
    for i=1:length(D)
        [mask, ~] = LMobjectmask(D(i).annotation,HOMEIMAGES);
        bw = logical(colorSegments(mask));
        bw = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));
        [~, name, ~] = fileparts(D(i).annotation.filename);
        imwrite(bw,strcat(HOMEMASKS,name,'_Binary.png'));
    end
    
    %LMdbshowobjects(LMquery(D, 'object.name','dn'),HOMEIMAGES);
    %LMdbshowobjects(LMquery(D, 'object.name','n'),HOMEIMAGES);
    
    % Extracting mask for Dim Nuclei (dn) Class labels
    [Ddn,~] = LMquery(D, 'object.name','nucleus');
    for i=1:length(Ddn)
        [mask, ~] = LMobjectmask(Ddn(i).annotation,HOMEIMAGES);
        bw = logical(colorSegments(mask));
        bw = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));
        [~, name, ~] = fileparts(D(i).annotation.filename);
        imwrite(bw,strcat(HOMEMASKS,name,'_Binary_Nucleus.png'));
    end
    
    [Ddn,~] = LMquery(D, 'object.name','dimnucleus');
    for i=1:length(Ddn)
        [mask, ~] = LMobjectmask(Ddn(i).annotation,HOMEIMAGES);
        bw = logical(colorSegments(mask));
        bw = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));
        [~, name, ~] = fileparts(D(i).annotation.filename);
        imwrite(bw,strcat(HOMEMASKS,name,'_Binary_DimNucleus.png'));
    end
    
    % Extracting mask for Normal Nuclei (n) class labels
    % query using 'n' is not working, might be single character comparison
    % problem. 
    % So I adopt different way to select normal nuclei (n)
    for i=1:length(Ddn)
        [mask, ~] = LMobjectmask(D(i).annotation,HOMEIMAGES);
        bw = logical(colorSegments(mask));
        bw1 = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));

        [mask, ~] = LMobjectmask(Ddn(i).annotation,HOMEIMAGES);
        bw = logical(colorSegments(mask));
        bw2 = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));
        
        bw = xor(bw1,bw2);
        [~, name, ~ ] = fileparts(D(i).annotation.filename);
        imwrite(bw,strcat(HOMEMASKS,name,'_Binary_DimNucleus.png'));
    end
end

