clear
clc

%refers to https://www.mathworks.com/help/vision/ug/object-detection-using-deep-learning.html#responsive_offcanvas

addpath(genpath('/Users/chunlidai/Documents/MATLAB/Examples/R2020b/'));
% addpath(genpath('/Users/chunlidai/surge/home/dai.56/arcticdemapp/landslide/code1/'));
addpath(genpath('/Users/chunlidai/ArchivedBoxSync/ESI2019/landslide/code1'));

% cifar10Data = tempdir;
% url = 'https://www.cs.toronto.edu/~kriz/cifar-10-matlab.tar.gz';
% helperCIFAR10Data.download(url,cifar10Data);
% [trainingImages,trainingLabels,testImages,testLabels] = helperCIFAR10Data.load(cifar10Data);
% size(testImages)

% % Prepare training images
% trainingLabels1=categorical({'slump', 'non slump'}');
% trainingLabels1(1)

            
numImageCategories = 2;

shapefiles='/Users/chunlidai/ArchivedBoxSync/AGU2019GSA2019/landslidecasesforTraining/Ashley/shapefiles/Slumps_2017_18_SCAR_FINAL.shp';
imagefiles={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};

[trainingImages1,trainingLabels1,slumpbox]=prepare(shapefiles,imagefiles);

categories(trainingLabels1)
% categories(testLabels)

[height,width,numChannels, ~] = size(trainingImages1);
imageSize = [height width numChannels];
inputLayer = imageInputLayer(imageSize);
filterSize = [5 5];
numFilters = 32;

middleLayers = [
% The first convolutional layer has a bank of 32 5x5x3 filters. A
% symmetric padding of 2 pixels is added to ensure that image borders
% are included in the processing. This is important to avoid
% information at the borders being washed away too early in the
% network.
convolution2dLayer(filterSize,numFilters,'Padding',2)
% Note that the third dimension of the filter can be omitted because it
% is automatically deduced based on the connectivity of the network. In
% this case because this layer follows the image layer, the third
% dimension must be 3 to match the number of channels in the input
% image.
% Next add the ReLU layer:
reluLayer()
% Follow it with a max pooling layer that has a 3x3 spatial pooling area
% and a stride of 2 pixels. This down-samples the data dimensions from
% 32x32 to 15x15.
maxPooling2dLayer(3,'Stride',2)
% Repeat the 3 core layers to complete the middle of the network.
convolution2dLayer(filterSize,numFilters,'Padding',2)
reluLayer()
maxPooling2dLayer(3, 'Stride',2)
convolution2dLayer(filterSize,2 * numFilters,'Padding',2)
reluLayer()
maxPooling2dLayer(3,'Stride',2)
];

finalLayers = [
% Add a fully connected layer with 64 output neurons. The output size of
% this layer will be an array with a length of 64.
fullyConnectedLayer(64)
% Add an ReLU non-linearity.
reluLayer
% Add the last fully connected layer. At this point, the network must
% produce 10 signals that can be used to measure whether the input image
% belongs to one category or another. This measurement is made using the
% subsequent loss layers.
fullyConnectedLayer(numImageCategories)
% Add the softmax loss layer and classification layer. The final layers use
% the output of the fully connected layer to compute the categorical
% probability distribution over the image classes. During the training
% process, all the network weights are tuned to minimize the loss over this
% categorical distribution.
softmaxLayer
classificationLayer
];

layers = [
    inputLayer
    middleLayers
    finalLayers
    ];
layers(2).Weights = 0.0001 * randn([filterSize numChannels numFilters]);

opts = trainingOptions('sgdm', ...
'Momentum', 0.9, ...
'InitialLearnRate', 0.001, ...
'LearnRateSchedule', 'piecewise', ...
'LearnRateDropFactor', 0.1, ...
'LearnRateDropPeriod', 8, ...
'L2Regularization', 0.004, ...
'MaxEpochs', 40, ...
'MiniBatchSize', 128, ...
'Verbose', true);

doTraining = true;
if doTraining
% Train a network.
cifar10Net = trainNetwork(trainingImages1, trainingLabels1, layers, opts);
else
% Load pre-trained detector for the example.
load('rcnnStopSigns.mat','cifar10Net')
end
% Extract the first convolutional layer weights
w = cifar10Net.Layers(2).Weights;
% rescale the weights to the range [0, 1] for better visualization
w = rescale(w);

figure
montage(w)

if 0 %classify a single image -> Yes or No output
% Run the network on the test set.
% testImages=trainingImages1(:,:,:,1);
idx=235-16:235+16-1;idy=1905-16:1905+16-1; % 235  1905; slump in slide 25
testImages=single(data0.z(idy,idx));

YTest = classify(cifar10Net, testImages);
% Calculate the accuracy.
% accuracy = sum(YTest == testLabels)/numel(testLabels);

end


% % Step 2 Box detection;

% data = load('stopSignsAndCars.mat', 'stopSignsAndCars');
% stopSignsAndCars = data.stopSignsAndCars;

slump1=slumpbox;

% Display one training image and the ground truth bounding boxes
I = imread(slump1.imageFilename{1});
I = insertObjectAnnotation(I,'Rectangle',slump1.slump{1},'slump','LineWidth',8);

% slump.slump(1)=[]; %x y dx dy
figure
imagesc(I)
colorbar;colormap jet;caxis([-10 10])

% A trained detector is loaded from disk to save time when running the
% example. Set this flag to true to train the detector.
doTraining = true;

if doTraining
    
    % Set training options
    options = trainingOptions('sgdm', ...
        'MiniBatchSize', 128, ...
        'InitialLearnRate', 1e-3, ...
        'LearnRateSchedule', 'piecewise', ...
        'LearnRateDropFactor', 0.1, ...
        'LearnRateDropPeriod', 100, ...
        'MaxEpochs', 100, ...
        'Verbose', true);
    
    % Train an R-CNN object detector. This will take several minutes.    
    rcnn = trainRCNNObjectDetector(slump1, cifar10Net, options, ...
    'NegativeOverlapRange', [0 0.3], 'PositiveOverlapRange',[0.5 1]);
else
    % Load pre-trained network for the example.
    load('rcnnStopSigns.mat','rcnn')       
end

% Read test image
changefile=imagefiles{1};
testImage = imread(changefile);

if  0 %Try it out on a test image

% Detect stop signs
% [bboxes,score,label] = detect(rcnn,testImage,'MiniBatchSize',128);
[bboxes,score,label] = detect(rcnn,testImage,'MiniBatchSize',32);

figure;imagesc(testImage);
colorbar;colormap jet;caxis([-10 10])
hold all
% ground truth
for i=1:28
clear box;box=slump1.slump(i);box=box{1}; %x y dx dy
x1=[box(1) box(1)+box(3)  box(1)+box(3)  box(1)  box(1)];
y1=[box(2) box(2)  box(2)+box(4)    box(2)+box(4)  box(2)];
plot(x1, y1, 'k.-','linewidth',3)
end
%classified results
for i=1:length(score)
clear box;box=bboxes(i,:); %x y dx dy
x1=[box(1) box(1)+box(3)  box(1)+box(3)  box(1)  box(1)];
y1=[box(2) box(2)  box(2)+box(4)    box(2)+box(4)  box(2)];
plot(x1, y1, 'r.-','linewidth',3)
end
end % if 0

% Get the classification scores map ^_^
featureMap = activations(rcnn.Network, testImage, 14); %layer 14 is the softmax layer.

% The softmax activations are stored in a 3-D array.
size(featureMap)
rcnn.ClassNames

stopSignMap = featureMap(:, :, 1); %Classification scores

[height, width, ~] = size(testImage);
stopSignMap = imresize(stopSignMap, [height, width]); %not accurate

% % Visualize the feature map superimposed on the test image. 
% featureMapOnImage = imfuse(testImage, stopSignMap); 
% 
% figure
% imshow(featureMapOnImage)

figure;imagesc(stopSignMap);colorbar;title('classification scores')
colorbar;colormap default;caxis([0 1])
hold all

M=stopSignMap>0.9;
B=bwboundaries(M,'noholes'); 
% plot the outlines
for k=1:length(B) %
    xid=B{k}(:,2); yid=B{k}(:,1);
% x=xout(xid);y=yout(yid);    
    xo{k}=xid;yo{k}=yid;
    hold on;plot(xid,yid,'r.-');
end

% figure;imagesc(testImage,'alphadata',M) %not showing well
% colorbar;colormap jet;caxis([-10 10])
% hold all


