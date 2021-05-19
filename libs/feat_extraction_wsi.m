function feat_extraction_wsi(folder_matpatches,folder_pyepistroma,folder_matcellmask,folder_savepath,quality)
placeholder = 'placeholder';
imgList=dir([folder_matpatches '*.mat']);
indx = randperm(numel(imgList));
numFiles=length(imgList);
imgList = {imgList(:).name};
qualitysaturationLim = quality.saturationLim;
qualityredChannelLim = quality.redChannelLim;
qualityblurLimit = quality.blurLimit;
qualityimgarea = quality.imgarea;
parfor nn=1:numFiles
    n= indx(nn);
    [~,imgName]=fileparts(imgList{n});
    % ToDo: Check if the nuclei and ES are already being processed, to be
    % used with another feature set
    imgFile=[folder_matpatches imgList{n}];
    outputFolder=[folder_savepath 'dataset_output/' imgName];
    maskFolder=[outputFolder '/png_binmask/png_cellmask/'];
    maskEFolder=[outputFolder '/png_binmask/png_cellepimask/'];
    maskSFolder=[outputFolder '/png_binmask/png_cellstromask/'];
    maskBFolder=[outputFolder '/png_binmask/png_cellbundmask/'];
    ESmaskFolder=[outputFolder '/png_binmask/png_epithmask/'];
    SSmaskFolder=[outputFolder '/png_binmask/png_stromask/'];
    EBmaskFolder=[outputFolder '/png_binmask/png_boundaryepistromask/'];
    featlocFolder=[outputFolder '/mat_simcellfeatures/'];
    disp(outputFolder)
    mkdir(outputFolder);
    mkdir(maskFolder);
    mkdir(maskEFolder);
    mkdir(maskSFolder);
    mkdir(maskBFolder);
    mkdir(ESmaskFolder);
    mkdir(SSmaskFolder);
    mkdir(EBmaskFolder);
    mkdir(featlocFolder);
    % open .mat file
    imgFile_load = load(imgFile);
    numtilestruct = length(imgFile_load.tileStruct);
    indx2 = randperm(numel(imgFile_load.tileStruct));
    fprintf('Processing file %s and tiles %d\n',imgName,numtilestruct);
    try
        curTile_ESmask_struct=h5read([folder_pyepistroma imgName '_ESmask'],'/mask');
    catch
        disp('disrupted or missing Epistroma Mask');
        filePh_error = fopen([folder_savepath 'errorlist_patches.txt'],'a');
        fprintf(filePh_error,'Disrupted epistroma mask %s with tiles %d\n',imgName,numtilestruct);
        fclose(filePh_error);
        % ToDo: delete the created folder
        rmdir(outputFolder,'s')
        continue;
    end
    curTile_Nmask_struct=h5read([folder_matcellmask imgName '_mask'],'/mask');
    for i=1:numtilestruct
        featFile=sprintf('%s/%s_%d.mat',featlocFolder,imgName,indx2(i));
        maskFile=sprintf('%s/%s_%d.png',maskFolder,imgName,indx2(i));
        maskEFile=sprintf('%s/%s_%d.png',maskEFolder,imgName,indx2(i));
        maskSFile=sprintf('%s/%s_%d.png',maskSFolder,imgName,indx2(i));
        maskBFile=sprintf('%s/%s_%d.png',maskBFolder,imgName,indx2(i));
        ESmaskFile=sprintf('%s/%s_%d.png',ESmaskFolder,imgName,indx2(i));
        SSmaskFile=sprintf('%s/%s_%d.png',SSmaskFolder,imgName,indx2(i));
        EBmaskFile=sprintf('%s/%s_%d.png',EBmaskFolder,imgName,indx2(i));
        curTile=imgFile_load.tileStruct(indx2(i)).data;
        if (getSaturationMetric(curTile)>qualitysaturationLim && ...
                getRedMetric(curTile)>qualityredChannelLim && ...
                blurMetric(curTile)>qualityblurLimit && ...
                getAreaTissue(curTile)>qualityimgarea)
            if exist(featFile,'file')~=2
                % Check if Cell mask and Epistroma mask have same quantity
                parsave(featFile, placeholder);
                if size(curTile_ESmask_struct) == size(curTile_Nmask_struct)
                    fprintf('Processing tile %s_%d\n',imgName,indx2(i));
                    curTile_ESmask = curTile_ESmask_struct(:,:,indx2(i))'; % For a single mask
                    curTile_Nmask = curTile_Nmask_struct(:,:,indx2(i))'; % For a single mask
                    [nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi, nucleiCentroids_bund,...
                        feat_simcell, feat_epi_simcell, feat_stro_simcell, feat_bund_simcell] = get_simpleCell(curTile,curTile_ESmask,curTile_Nmask,maskFile,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile);
                    parsave_cellfeat(featFile,nucleiCentroids,nucleiCentroids_stro,nucleiCentroids_epi,nucleiCentroids_bund,...
                        feat_simcell,feat_epi_simcell,feat_stro_simcell,feat_bund_simcell);
                else
                    filePh_error = fopen([folder_savepath 'errorlist_patches.txt'],'a');
                    fprintf(filePh_error,'Mismatch cell and epistroma mask %s_%d\n',imgName,indx2(i));
                    fclose(filePh_error);
                    continue
                end
            else
                %disp('mat tile already exists')
                continue
            end
        end
    end
end
%ToDo: include method to delete the 1Kb files
parfor n=1:numFiles
    [~,imgName]=fileparts(imgList{n});
    % ToDo: Check if the nuclei and ES are already being processed, to be
    % used with another feature set
    %imgFile=[folder_matpatches imgList{n}];
    outputFolder=[folder_savepath 'dataset_output/' imgName];
    featlocFolder=[outputFolder '/mat_simcellfeatures/*.mat'];
    files = dir(featlocFolder);
    for ii = 1:length(files)
        if files(ii).bytes==183 % file with 'placeholder' word
            delete(fullfile(files(ii).folder, files(ii).name))
        end
    end
end
fprintf('Done!\n');
% Save the geenral feature description
numZernikePol=36;
%zernNames={};
zernNames = cell(1,numZernikePol*2);
zernumb = 0;
for i=1:2:numZernikePol*2
    zernumb = zernumb+1;
    zernNames{i}= ['ZernPol' num2str(zernumb) '_A'];
    zernNames{i+1}= ['ZernPol' num2str(zernumb) '_Phi'];
end
harNames={'AngularSecondMoment','Contrast','Correlation',...
    'Variance','Homogeneity','SumAverage','SumVariance','SumEntropy',...
    'Entropy','DifferenceVariance','DifferenceEntropy','InfoMeasureCorrI',...
    'InfoMeasureCorrII','MaxCorrCoeff'};
description_simcell=[{'Area','Eccentricity','RatioAxes','MedianRed','EntropyRed',...
    'MinIntensity','MaxIntensity','EquivDiameter','Orientation',...
    'EntropyIntensity','MedianIntensity','EdgeMeanIntensity','RatioMedianRedBlue',...
    'RatioMedianRedGreen'},harNames,zernNames];
save([folder_savepath 'simcell_feature_description.mat'],'description_simcell');
end