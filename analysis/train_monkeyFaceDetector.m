%% train_monkeyFaceDetector
%

clc;
clear all;

%% setup
% path
imDir = 'B:\nicolas.borderies\Projets scientifiques\Projet PNH\Pharmaco-Motiv-actual\images';
addpath(imDir);

% load positive image sample 
positiveFolder = [imDir '\macaque face database'];
cd(positiveFolder);
load('face_data.mat');

% load negative image sample
negativeFolder = [imDir '\negative_image_google'];

%% training

trainCascadeObjectDetector('monkeyFaceDetector.xml',face_data,negativeFolder,'FalseAlarmRate',0.2,'NumCascadeStages',5);

%% testing

detector = vision.CascadeObjectDetector('monkeyFaceDetector.xml');

% Read the test image.
img = imread('img_test2.jpg');

% Detection
bbox = step(detector,img);

% Insert bounding boxes and return marked image.
detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');

% Display the detected stop sign.
figure;
imshow(detectedImg);