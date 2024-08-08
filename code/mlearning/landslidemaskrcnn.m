% He et al., 2017
% example from https://github.com/matlab-deep-learning/mask-rcnn
% see /Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/mask-rcnn-main/MaskRCNNTrainingExample.mlx

clear
clc

addpath(genpath('/Users/chunlidai/Documents/MATLAB/Examples/R2020b/'));
addpath(genpath('/Users/chunlidai/surge/home/dai.56/arcticdemapp/landslide/code1/'));
addpath(genpath('/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/'));

% Set directory locations for the coco dataset
trainImgFolder = '/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/mask-rcnn-main/coco/train2014';
annotationfile = '/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/mask-rcnn-main/coco/annotations/instances_train2014.json';


% Define teh categories to train the network on
trainCats = {'person','car'};

unpackAnnotationFolder = '/Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/mask-rcnn-main/coco/annotations_unpacked/matFiles';

helper.unpackAnnotations(trainCats, annotationfile, trainImgFolder, unpackAnnotationFolder);

% Set image size for training
imageSize = [800 800 3];

% Create the training datastore to read image and ground truth data from
% the unpacked annotation MAT file
ds = fileDatastore(unpackAnnotationFolder, 'ReadFcn',@(x)helper.cocoAnnotationMATReader(x, trainImgFolder));

trainDS = transform(ds, @(x)helper.preprocessData(x, imageSize));

trainDS.shuffle();
data = preview(trainDS);

classNames = trainCats;
numClasses = numel(classNames);
% Add a background class
classNames = [classNames {'background'}];

params = createMaskRCNNConfig(imageSize, numClasses, classNames);
disp(params);

if canUseGPU
    executionEnvironment = "gpu";
else
    executionEnvironment = "cpu";
end

% Download and load the pre-trained network -

doTraining = false;
if ~doTraining
    datadir = tempdir;
    url = 'https://www.mathworks.com/supportfiles/vision/data/maskrcnn_pretrained_person_car.mat';
    
    helper.downloadTrainedMaskRCNN(url, datadir);
    pretrained = load(fullfile(datadir, 'maskrcnn_pretrained_person_car.mat'));
    net = pretrained.net;
    
    % Extract Mask segmentation sub-network
    maskSubnet = helper.extractMaskNetwork(net);
else

% Create the maskRCNN network
dlnet = createMaskRCNN(numClasses, params);

% SGDM learning parameters
initialLearnRate = 0.01;
momemtum = 0.9;
decay = 0.0001;
velocity = [];
maxEpochs = 30;

minibatchSize = 2;

myMiniBatchFcn = @(img, boxes, labels, masks) deal(cat(4, img{:}), boxes, labels, masks);

mb = minibatchqueue(trainDS, 4, "MiniBatchFormat", ["SSCB", "", "", ""],...
                            "MiniBatchSize", minibatchSize,...
                            "OutputCast", ["single","","",""],...
                            "OutputAsDlArray", [true, false, false, false],...
                            "MiniBatchFcn", myMiniBatchFcn,...
                            "OutputEnvironment", [executionEnvironment,"cpu","cpu","cpu"]);
                        
numEpoch = 1;
numIteration = 1; 

start = tic;
if doTraining
    
     % Create subplots for the learning rate and mini-batch loss.
    fig = figure;
    [lossPlotter] = helper.configureTrainingProgressPlotter(fig);
    
    % Initialize verbose output
    helper.initializeVerboseOutput([]);
    
    % Custom training loop.
    while numEpoch < maxEpochs
    mb.reset();
    mb.shuffle();
    
        while mb.hasdata()
            % get next batch from minibatchqueue
            [X, gtBox, gtClass, gtMask] = mb.next();
        
            % Evaluate the model gradients and loss using dlfeval
            [gradients, loss, state] = dlfeval(@networkGradients, X, gtBox, gtClass, gtMask, dlnet, params);
            dlnet.State = state;
            
            % compute the learning rate for current iteration
            learnRate = initialLearnRate/(1 + decay*numIteration);
            
            if(~isempty(gradients) && ~isempty(loss))
    
                [dlnet.Learnables, velocity] = sgdmupdate(dlnet.Learnables, gradients, velocity, learnRate, momemtum);
            else
                continue;
            end
            helper.displayVerboseOutputEveryEpoch(start,learnRate,numEpoch,numIteration,loss);
                
            % Plot loss/ accuracy metric
             D = duration(0,0,toc(start),'Format','hh:mm:ss');
            addpoints(lossPlotter,numIteration,double(gather(extractdata(loss))))
            
            subplot(2,1,2)

            title("Epoch: " + numEpoch + ", Elapsed: " + string(D))
           
            drawnow
            
            numIteration = numIteration + 1;
    
        end
    numEpoch = numEpoch + 1;
    
    end
end

net = dlnet;
maskSubnet = helper.extractMaskNetwork(net);
end %if training

% Read the image for inference
img = imread('visionteam.jpg');

% Define the target size of the image for inference
targetSize = [700 700 3];

% Resize the image maintaining the aspect ratio and scaling the largest
% dimension to the target size.
imgSize = size(img);
[~, maxDim] = max(imgSize);
resizeSize = [NaN NaN]; 
resizeSize(maxDim) = targetSize(maxDim);

img = imresize(img, resizeSize);

% detect the objects and their masks
[boxes, scores, labels, masks] = detectMaskRCNN(net, maskSubnet, img, params, executionEnvironment);

% Visualize Predictions

% Overlay the detected masks on the image using the insertObjectMask
% function.
if(isempty(masks))
    overlayedImage = img;
else
    overlayedImage = insertObjectMask(img, masks);
end
figure, imshow(overlayedImage)

% Show the bounding boxes and labels on the objects
showShape("rectangle", gather(boxes), "Label", labels, "LineColor",'r')


%%

% Perform prediction using the detectMaskRCNN function-

if 0 %use the trained net
    datadir = tempdir;
    url = 'https://www.mathworks.com/supportfiles/vision/data/maskrcnn_pretrained_person_car.mat';
    
    helper.downloadTrainedMaskRCNN(url, datadir);
    pretrained = load(fullfile(datadir, 'maskrcnn_pretrained_person_car.mat'));
    net = pretrained.net;
    
    % Extract Mask segmentation sub-network
    maskSubnet = helper.extractMaskNetwork(net);    
    
img = imread('visionteam.jpg');
executionEnvironment = "cpu";

% helper.createNetworkConfiguration %to check
[height,width,numChannels, ~] = size(img);
imageSize = [height width numChannels];
numClasses=2;
classNames={'Foreground', 'Background'};
params = createMaskRCNNConfig(imageSize, numClasses, classNames);

[boxes, scores, labels, masks] = detectMaskRCNN(net, maskSubnet, img, params, executionEnvironment);

% Visualize results
overlayedImage = insertObjectMask(img, masks);
figure, imshow(overlayedImage)
showShape("rectangle", gather(boxes), "Label", labels, "LineColor",'r')

end
% To train your Mask-RCNN network, follow the steps outlined in the following examples.

% MaskRCNNTrainingExample.mlx 

