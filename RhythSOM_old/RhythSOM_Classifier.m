function [clusData, varargout] = RhythSOM_Classifier( Data, minPCvar, autoClus, numRep, waveletParams, dirSave)

%    RhythSOM_Classifier is RhythSOM main function. It is thought to display
% clusterized electrophysiological rhythms, such hippocampal sharp wave
% ripples or theta cycles, by their waveform. A 'Data' matrix with each
% sample in a row first goes through a Principal Components Analysis (PCA), 
% and the selected PCs are fed into a Self-Organizing Map (SOM) algorithm from
% the SOM Toolbox (link at the bottom of this description). Units from the
% map are clusterized by a k-means algorithm, and results of these 
% clusterization are finally plotted in two separate figures. On the left,
% for each cluster mean waveforms and their respective wavelet, and on the
% right, the Self-Organized Map with the mean wavelet on top. 
%
%    There are two parameters that can set different options:
% 1) The number of Principal Components can be chosen manually ('minPCvar=nan'),
%    or automatically by selecting the minimum explained variance with 'minPCvar'
% 2) The number of clusters can be chosen also manually ('autoClus=0'), or
%    automatically ('autoClus=1').
% 
% 
% clusData, [clussMap, sMap] = RhythmSOM_Classifier( Data, [minPCvar, autoClus, numRep, waveletParams, imageDir] )
% 
% Inputs
%
%   Data          (matrix)  N x M matrix, formed by N samples of size M. 
%                           Examples of usable data could be N 100ms-long
%                           somatic sharp wave ripples, or N theta cycles.
%
%   minPCvar      (num/nan) Optional. Minimum explained variance for the
%                           selected number of Principal Components. For
%                           example, if 'minPCvar=0.99',the number of PCs
%                           will be chosen automatically to be the one that
%                           explains 99% of variance. If 'minPCvar=nan',
%                           then the number of PCs will be chosen manually
%                           based on a graph. Variable can be then 'nan' or
%                           a float between 0 and 1.
%
%   autoClus        (bool)  Optional. Select number of clusters automatically
%                           (1) or manually (0)? Automatically by default.
%
%   numRep           (int)  Optional. Number of repetitions for k-means
%                           algorithm. If numRep is greater than 1, several
%                           iterations will be made, and results from best
%                           iteration will be selected. This way the performance
%                           doesn't depend so much on local minima near the
%                           first random centroids.
%
%   waveletParams (struct)  Optional. Structure with parameters to plot the 
%                           wavelet on Figure 1: sample frecuency first, 
%                           frequency limits second, and number of consecutive
%                           waveforms to be displayed. If not given,
%                           wavelets are not shown.
%                           Eg.: 20000 Hz, frequencies between 15 and 250 Hz
%                                and two repeated theta cycles displayed
%
%                                   waveletParams = {20000, [15 250], 2}
%
%                           Eg.: 5000 Hz, frequencies between 70 and 500 Hz
%                                and ripple is displayed just once
%  
%                                   waveletParams = {5000, [70 500], 1}
%
%   dirSave       (string)  Optional. Directory to save the ouput figures.
%                           If not given, they are saved on the current
%                           location.
%                           Eg.: '/home/Projects/RhythSOM/images/'
% 
% 
% 
% Outputs
%
%   clusData     (matrix)  1 x N matrix, N being the number of samples.
%                          For each sample, the cluster number to which
%                          it belongs. 
%                          Eg.: clusData = 
%                                   1   1   1   2   2   3   3   3   3
%
%   clussMap     (matrix)  Optional. 1 x U matrix, U being the number of units
%                          on the sMap. For each unit, the cluster number to
%                          which it belongs.
%                          Eg.: clussMap = 
%                                   1   1   2   3   3
%
%   sMap         (struct)  Optional. SOM ouput sctructure from som_make.
%
%   Figure 1     (figure)  'RhythSOM_meanWaveforms.png'. Save on dirSave if
%                          given; if not, in the current folder.
%
%   Figure 2     (figure)  'RhythSOM_map.png'. Save on dirSave if given;
%                          if not, in the current folder.
%
%
%
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%
%                      
% LCN-acnavasolive 2019

    % 'saveName' optional variable
    if nargin < 6; dirSave = ''; end
    % 'waveletParams' optional variable
    if nargin < 5; waveletParams = {}; end
    % 'autoClus' optional variable
    if nargin < 4; numRep = 100; end 
    % 'autoClus' optional variable
    if nargin < 3; autoClus = 1; end 
    % 'minPCvar' optional variable
    if nargin < 2; minPCvar = 0.8; end
        

    % 1. Principal Components matrix of Data. 
    %    Manual or automatic selection of number of PCs depending on 'minPCvar'
    DataPCA = RhythSOM_pcs(Data, minPCvar);
    
    % 2. Normalized DataPCA so SOM algorithm works with normalized 
    %    euclidean distances, and avoid biases.
    [DataNorm, DataDenorm] = RhythSOM_normdata(DataPCA);

    % 3. Training Self-Organized-Map
    sMap = som_make(DataNorm,'lattice','rect','training','long','tracking',1);
    % Find the best-matching units (bmus) from the map for the given vectors
    [BMUs , BMUerror] = som_bmus(sMap, DataNorm);
    
    % 4. Clusterization by k-means
    [clusData, clussMap, clusNum] = RhythSOM_clusters(sMap, BMUs, autoClus, numRep);
    
    % 5. Visualization of clustering sMap
    RhythSOM_plot(Data, sMap, clusData, clussMap, waveletParams, dirSave)
    
    % Optional ouputs
    %  - Clusters of SOM
    varargout{1} = clussMap;
    %  - Self organized map
    varargout{2} = sMap;
    
end