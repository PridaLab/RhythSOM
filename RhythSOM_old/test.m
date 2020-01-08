% Script to test RhythSOM_Classifier
clear all; clc; close all;

% Add all scripts to our path
addpath(genpath(pwd));

% Test data
Data = [ 2+rand(200,20), 4+rand(200,20), 0+rand(200,20), 0+rand(200,20); ...
         0+rand(200,20), 0+rand(200,20), 2+rand(200,20), 3+rand(200,20)]';
     
% Define parameters
%  - Minimum explained variace for Principal Components Analysis
minPCvar = 0.9;
%  - Automatic clusterization (1: true, 0: false)
autoClus = 1;
%  - Number of repetitions of the k-means algorithm to clusterize 
numRep = 50;
%  - Parameters to plot wavelet: sampling frequency, frequency limits and
%    number of consecutive waveforms displayed
waveletParams = {20000, [60 250], 1};
%  - Directory to save figures
dirSave = '/home/andrea/Projects/RhythSOM/images/';

% Ripple matrix
rippleMat = load('/home/andrea/Projects/RhythSOM/data/db_mat_ripples.mat');
rippleMat = rippleMat.db_mat_ripples;
rippleExp = fieldnames(rippleMat);
Data = [];
for ii = 1:length(rippleExp)
    if contains(rippleExp{ii}, '5verde')
        idxs = (20001-1)/2 - 0.045*waveletParams{1} : (20001-1)/2 + 0.045*waveletParams{1};
        Data = [Data; rippleMat.(rippleExp{ii}).spont_mat_ripple(idxs,:)'];
        idxs = 5100 : 5100 + 0.090*waveletParams{1};
        Data = [Data; rippleMat.(rippleExp{ii}).evo_mat_ripple(idxs,:)'];
    end
end

% Run
close all
clusData = RhythSOM_Classifier( Data, minPCvar, autoClus, numRep, waveletParams, dirSave);













