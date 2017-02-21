function ComputeFeaturesSummary(FeaturePath)
%
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    srcFiles1 = dir(strcat(FeaturePath,'*_CMFeatures.csv'));
	srcFiles2 = dir(strcat(FeaturePath,'*_RLFeatures.csv'));
	srcFiles3 = dir(strcat(FeaturePath,'*_IntensityFeatures.csv'));
	srcFiles4 = dir(strcat(FeaturePath,'*_MorphologicalFeatures.csv'));

	SumFeatures = table;
	SumFeatures_2 = table;
	ImageNames = cell(length(srcFiles1),1);
	for i = 1:length(srcFiles1)

		[~, ImageName, ~] = fileparts(srcFiles1(i).name);
		ImageNames{i} =  ImageName(1:end-11);
		
		 CM = readtable(strcat(FeaturePath,srcFiles1(i).name));
		 %CM = CM(:,1:14);
		 [~,NCol] = size(CM);
		 CM = CM(:,1:NCol-1);
		 RL = readtable(strcat(FeaturePath,srcFiles2(i).name));
		 %RL = RL(:,1:11);
		 [~,NCol] = size(RL);
		 RL = RL(:,1:NCol-1);
		 In = readtable(strcat(FeaturePath,srcFiles3(i).name));
		 %In = In(:,1:9);
		 [~,NCol] = size(In);
		 In = In(:,1:NCol-1);
		 Mo = readtable(strcat(FeaturePath,srcFiles4(i).name));
		 Mo = Mo(:,1:13);
	 
		csize=size(CM);
		msize=size(Mo);
		less = csize > msize;
		if less(1)
			CM = CM(1:msize,:);
			RL = RL(1:msize,:);
		end
		less = csize < msize;
		if less(1)
			Mo = Mo(1:csize,:);
			In = In(1:csize,:);
		end
		
		Features = horzcat(Mo,In,CM,RL);
		writetable(Features,strcat(FeaturePath,ImageNames{i},'_AllFeatures.csv'));            
		
		SumFeatures(i,:) = SummarizeFeatures(Features);
		SumFeatures_2(i,:) = SummarizeFeatures_2(Features);
	end
	SumFeatures = horzcat(cell2table(ImageNames),SumFeatures);
	writetable(SumFeatures,strcat(FeaturePath,'Features_Exp.csv'));

	SumFeatures_2 = horzcat(cell2table(ImageNames),SumFeatures_2);
	writetable(SumFeatures_2,strcat(FeaturePath,'Features_Exp_Summary2.csv'));
end