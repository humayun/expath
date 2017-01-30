function GenerateMaskFromAnnotation()

    %% Downloadiing User Database (Images + Annotations) from LabelMe Website
    % Define the root folder for the images
    HOMEIMAGES = 'H:/Datasets/LSM/LabelMe/PreExpanded'; 
    HOMEANNOTATIONS = 'H:/Datasets/LSM/LabelMe/PreExpanded/Annotation'; 
    HOMEMASKS = 'H:/Datasets/LSM/LabelMe/PreExpanded/Masks/'; 
    mkdir(HOMEMASKS);

    % Database name with path on LabelMe website
    DatabaseName ={'users/hirshad/becklab_perexpanded'};

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

for i=1:length(Files)
    [~,ImageName,~] = fileparts(Files(i).name);
    im = imread(strcat(BWPath, ImageName, '.png'));
    ImageName(ImageName==' ')='';
    imwrite(im, lower(strcat(BWPath2, ImageName, '.png')));
end
