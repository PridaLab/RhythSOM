function [DataNorm, DataDenorm] = RhythSOM_normdata(DataPCA)


%   RhythSOM_normdata computes normalized DataPCA so SOM algorithm works 
% with normalized euclidean distances, and avoid biases.
%
% [DataNorm, DataDenorm] = RhythSOM_normdata(DataPCA)
%
% Inputs
%
%   DataPCA  (matrix)  N x K matrix, formed by N samples of size K 
%                      (selected number of PCs). It's the Principal
%                      Components matrix of Data.
%
% Ouputs
%
%   DataNorm (struct)  Structured data to feed som_make().
%                      DataNorm.data is a N x K normalized DataPCA matrix
%                      so SOM algorithm can work with normalized euclidean
%                      distances, and avoid biases.
%
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%
%                      
% LCN-acnavasolive 2019

    % Structure: data, 'name'+name, 'comp_names'+component names
    DataPCA = som_data_struct(DataPCA,'name','HippoRhythms');

    % SOM is based on Euclidian distances, we normalize components
    DataNorm = som_normalize(DataPCA,'var');
    x = DataNorm.data(1,:);

    % Denormalization: Denormalized data
    DataDenorm = som_denormalize(x,DataNorm);
 
end