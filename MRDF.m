%%% This is the code for MRDF produced by Junxiong Zhou.
%%% Email: zhou1743@umn.edu
%%% Copyright belongs to Junxiong Zhou
%%% When using the code, please cite the following paper
%%% Zhou, J., Qiu, Y., Chen, J., & Chen, X. (2021). A geometric misregistration resistant data fusion approach for adding red-edge (RE) and short-wave infrared (SWIR) bands to high spatial resolution imagery. Science of Remote Sensing, 100033.
%%% Matlab >= R2020a, https://www.mathworks.com/help/map/ref/readgeoraster.html;
%%% Feel free to contact me if you have any question.

clc;clear;close;
addpath('.\functions\');

%% Before using the code, please set the following parameters according to your dataset
path = 'test data\'; % data path

VNIR = [1,2,3,7];  % VNIR bands of the medium image
SWIR = [8,9];     % SWIR bands of the medium image

coarse = double(readgeoraster([path 'Medium'])); % Open a Medium resolution image
fine = double(readgeoraster([path 'Fine']));     % Open a Fine resolution image
OutName = [path 'MRDF.tiff'];            % output file name

mw = 5;          % size of moving window for searching similar pixels
ExR = 2;         % upscaling parameter for extended box
numM = 10;       % number of similar pixels

[nsc,nlc,nbc] = size(coarse);
[ns,nl,nb] = size(fine);
ratio = ns/nsc*1.0; % ns/nsc should be equal to nl/nlc

%% local model: similar pixel searching
mw = mw*ratio;
PreFine = zeros(ns,nl,nbc)*1.0;
parfor i = 1:ns
    for j = 1:nl
        ai = max(1, i-mw); bi = min(ns, i+mw);
        aj = max(1, j-mw); bj = min(nl, j+mw);

        aic = floor(ai/ratio)+1; bic = floor(bi/ratio);
        ajc = floor(aj/ratio)+1; bjc = floor(bj/ratio);

        NP = reshape(coarse(aic:bic,ajc:bjc,:),(bic-aic+1)*(bjc-ajc+1),nbc);  
        candidateP = NP(:,VNIR);     

        center = fine(i,j,:); center = center(:);
        EuDis = sqrt(sum((center - candidateP').^2,1))';

        num = min(length(EuDis), numM);
        [EuDisS,indB] = sort(EuDis);
        indB = indB(1:num);

        weight = EuDisS(1:num)+1e-5; weight = 1./weight/sum(1./weight);
        similarP = (NP(indB,:)'*weight)';
        coefficient=[similarP(1,VNIR);ones(1,nb)]'\center;
        index = find(isnan(coefficient));

        F = similarP(1,:)*coefficient(1)+coefficient(2);
        PreFine(i,j,:) = F;
    end
end

%% global model: partial linear regression
w=1; sigma=ratio/2;
PSF=PSF_template(ratio,w,sigma);
PAN_upscaled=dowmsample_cube(fine,ratio,w,PSF);
GB0=D3_D2(PAN_upscaled);
GB1 = D3_D2(coarse);
[~,~,~,~,xrc1,~] = plsregress(GB0',GB1',nb-1);
GBF=D3_D2(fine);
Ff1=[ones(ns*nl,1),GBF']*xrc1; 
PreFine_global = reshape(Ff1, ns, nl, nbc);

clearvars PAN_upscaled GB0 GB1 GBF Ff1 %clear temporary variables

%% combing local and global predictions for SWIR bands
localSWIR = PreFine(:,:,SWIR);
globalSWIR = PreFine_global(:,:,SWIR);

s = ExR*ratio; PSFh=PSF_template(ExR,1,1);
coarseSWIR = dowmsample_cube(coarse(:,:,SWIR),ExR,w,PSFh);
sigma=s/2; PSFh=PSF_template(s,w,sigma);
RB=coarseSWIR-dowmsample_cube(localSWIR,s,w,PSFh);
RB_global = coarseSWIR - dowmsample_cube(globalSWIR,s,w,PSFh);

residualsC = zeros(nsc/ExR,nlc/ExR,length(SWIR));
PreSWIR = zeros(ns,nl,length(SWIR));
for band = 1:length(SWIR)
    weight = (1./abs(RB(:,:,band))+1e-8) ./ (1./abs(RB(:,:,band))+1./abs(RB_global(:,:,band))+1e-8);
    weightR = imresize(weight, s, 'nearest');
    residualsC(:,:,band) = weight.*RB(:,:,band) + (1-weight).*RB_global(:,:,band);
    PreSWIR(:,:,band) = weightR.*localSWIR(:,:,band) + (1-weightR).*globalSWIR(:,:,band);
end

residualsC = BlockEffects(residualsC, dowmsample_cube(fine,s,w,PSFh), 1, 4); %spatial filtering
residuals = imresize(residualsC, s, 'bicubic');

PreSWIR = PreSWIR + residuals;
PreSWIR = BlockEffects(PreSWIR, fine, 5, min(0.5*ratio,30)); %spatial filtering
PreFine(:,:,SWIR) = PreSWIR;

clearvars localSWIR globalSWIR residualsC residuals PreSWIR RB RB_global PreSWIR %clear temporary variables

%% output
t = Tiff(OutName,'w');
tagstruct.ImageLength = ns; 
tagstruct.ImageWidth = nl;  
tagstruct.Photometric = 1;
tagstruct.BitsPerSample = 64;
tagstruct.SamplesPerPixel = nbc;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB'; 
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
t.setTag(tagstruct);

t.write(PreFine);
t.close;




