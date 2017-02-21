function [ p ] = InitializeParameters( )
%InitializeParameters inilialize parameters function
%   It inilialize parameters of nuclei segmentation and feature extraction 
% pipeline
%
% Author: (12/2015)
% -------------------------------------------
% Humayun Irshad (humayun.irshad@gmail.com)
% BIDMC, Harvard Medical School
% -------------------------------------------

    % General Parameters 
    p.spacing = [1,1];
    p.debugMode = false;
    p.poolSize = 1;
    p.flagParallelize = false;

    % Foreground Segmentation Parameters
    p.threshmodel = 'poisson';
    p.localWindowRadius = 100;
    p.localWindowPace = p.localWindowRadius/3;
    p.minLocalGlobalThresholdRatio = 0.60;
    %low value of minLocalGlobalThresholdRatio generate more foreground regions
    p.numHistogramBins = 256;

    % Seed Detection Parameters
    p.maxNumLoGScales = 10; 
    p.flagBrightBlobs = true;    % false
    p.logResponseCutoff = eps;   % 0
    p.showPlots = false;         
    p.flagShowNucleiSeedMask = true;
    p.flagShowNucleiSizeMask = false;

    % Nuceli Selection Parameters
    p.thresholdvalue = 0.001;
    p.minNucleiArea = 10000;
    p.maxNucleiArea = 100000;
    p.nucleiDiameterRange = [125,275];
    p.meanNucleiDiamter = mean(p.nucleiDiameterRange);
    p.meanNucleiRadius = 0.5 * p.meanNucleiDiamter;
    p.minDistanceBWSeeds = 0.5 * p.meanNucleiDiamter;
    p.minNucleiBBox = 0.3 * min(p.nucleiDiameterRange);
    p.seedToBackgroundDistanceCutoff = 0.2 * min(p.nucleiDiameterRange);
    %low value of seedToBackgroundDistanceCutoff(0.3) generate more seed, 
    %while high value (0.5) generate less seed points

    % Features Computation
    p.GrayLevels = 256;
end

