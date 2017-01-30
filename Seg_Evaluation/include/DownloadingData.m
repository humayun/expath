clear all;

%% Downloadiing User Database (Images + Annotations) from LabelMe Website
% Define the root folder for the images
HOMEIMAGES = 'H:/Datasets/LSM/LabelMe/Matlab/Images'; 
HOMEANNOTATIONS = 'H:/Datasets/LSM/LabelMe/Matlab/Annotation'; 

% Database name with path on LabelMe website
DatabaseName ={'users/hirshad/test'};

LMinstall(DatabaseName, HOMEIMAGES,HOMEANNOTATIONS);

%% Visualization of Annotation
% Reading the index
D = LMdatabase(HOMEANNOTATIONS);

% Displaying the images/objects
LMplot(D,1,HOMEIMAGES);
LMdbshowscenes(D(1:4),HOMEIMAGES);
LMdbshowobjects(D(1),HOMEIMAGES);

% Displaying the object statistics in database
LMobjectnames(D);
[names, counts] = LMobjectnames(D);

%% Extracting polygon
[x,y] = LMobjectpolygon(D(1).annotation,'n');

figure
plot(x{1},y{1},'r')
axis('ij')

%% Extracting segments
[mask, class] = LMobjectmask(D(1).annotation,HOMEIMAGES);
figure, imshow(colorSegments(mask));

bw = logical(colorSegments(mask));
bw = max(max(bw(:,:,1),bw(:,:,2)), bw(:,:,3));
figure, imshow(bw);
