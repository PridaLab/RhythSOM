function [clusData, clussMap, varargout] = RhythSOM_probabilisticClusters(sMap, BMUs, method, numTrials)


%   RhythSOM_probabilisticClusters clusterizes the sMap Best Matching
% Units (BMUs) through a variation of the k-means method. The function
% uses RhythSOM_clusters over numTrial times to create an association
% probability matrix between sMap's units and clusterize themin the most
% probable way. 
% 
%   [clusData, clussMap, clusNum] = RhythSOM_probabilisticClusters(sMap, BMUs, method, numTrials)
% 
% Inputs
%
%   sMap     (struct)  SOM ouput sctructure from som_make.
%
%   BMUs     (matrix)  Best-matching units output matrix from som_bmus
%                      for a given sMap.
%
%   method   (string)  Kmeans method used, 'batch' or 'seq'.
%
%   autoClus   (bool)  Optional. Select number of clusters automatically
%                      (1) or manually (0)? Automatically by default.
%
%   numTrials   (int)  Optional. Number of repetitions for 
%                      RhythSOM_clusters execution to construct
%                      the association probability matrix. The greater 
%                      the matrix, the greater the fiability. A minimum
%                      of 100 is recommended.
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
%   clusterMap   (matrix)  1 x U matrix, U being the number of units on 
%                          the sMap. For each unit, the cluster number to
%                          which it belongs.
%                          Eg.: clussMap = 
%                                   1   1   2   3   3
%
%   clusNum         (int)  Optional. Number of clusters.
%
%   W            (matrix)  Optional. U x U matrix, U being the number of
%                          units on the sMap. The values represent the
%                          probability of association between two units.
%
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%                          
% LCN - E.F Rodríguez Sebastián 2019

    % 'method' optional variable
    if nargin < 3; method = 'batch'; end
    % 'numTrials' optional variable
    if nargin < 4; numTrials = 100; end 
    
    %Results matrix creation
    trialClusters = [];
    trialNumClusters = [];
    trialClusMap = [];
    trialClusCentroids = [];

    %Test execution
    for trial = 1:numTrials
        [clusData, clussMap, clusNum, clusCentroids] = RhythSOM_clusters(sMap, BMUs, method, 1, 100);
        close
        trialClusters(:,trial) = clusData;
        trialNumClusters(1,trial) = clusNum;
        trialClusMap(:,trial) = clussMap;
        try
            trialClusCentroids(:,:,trial) = clusCentroids;
        catch
            trialClusCentroids(:,:,trial) = NaN(size(trialClusCentroids(:,:,trial-1)));
        end
        clc
        fprintf('Finished trial (%d/%d)\n', trial, numTrials)
    end
    clear clusNum

    %Look for cases with different numbers of clusters and take the majoritary case
    cases = unique(trialNumClusters);
    casesHits = [];
    for c = 1:length(cases)
        casesHits(c) = length(find(trialNumClusters == cases(c)));
    end
    clusNum = cases(find(casesHits == max(casesHits)));
    delInd = find(trialNumClusters ~= clusNum);
    trialClusters(:,delInd) = [];
    trialClusMap(:,delInd) = [];
    trialClusCentroids(:,:,delInd) = [];

    %Computation of association matrix (number of times each unit is grouped with other units)
    W = [];
    for unit = 1:size(trialClusMap, 1)
        vec = zeros(1,size(trialClusMap, 1));
        for trial = 1:size(trialClusMap, 2)
            clus = trialClusMap(unit, trial);
            coincidences = find(trialClusMap(:, trial) == clus);
            vec(coincidences) = vec(coincidences) + 1;
        end
        W(unit, :) = vec;
    end
    W = W / size(trialClusMap, 2);

    %Look for templates of the N (clusNum) units with most associations
    A = W;
    A(A<1) = 0;
    U = unique(A, 'rows');
    sumU = sum(U, 2);
    [~,I] = sort(sumU);
    idx = I((end-(clusNum-1)):end);
    clusTemplates = U(idx, :);

    %Create the cluster map vector and give values to the units that match the
    %corresponding template
    clussMap = NaN(size(W,1), 1);
    for cTmp = 1:size(clusTemplates, 1)
        for unit = 1:size(A, 1)
            if sum(abs(clusTemplates(cTmp,:) - A(unit,:))) == 0
                clussMap(unit) = cTmp;
            end
        end
    end

    %Append the values that are missing using the greatest probability of 
    %its row
    nonClusSet = find(isnan(clussMap));
    N = length(nonClusSet) * 4;
    P = W;
    for unit = 1:size(P, 1)
        P(unit,unit) = 0;
        if ~ismember(unit, nonClusSet)
            idx = find(P(unit,:) == 1);
            P(unit,idx) = 0;
        end
    end
    %Check association probabilities of each unit in the non clustered
    %set. The units will be appended to the cluster of the max association
    %probability. If after a stop criteria (based on the initial number
    %of elements to clusterize) there are unclustered units, a residual 
    %cluster will be created.
    i = 0;
    checkVec = [];
    while sum(isnan(clussMap)) ~= 0
        if i < length(nonClusSet)
            i = i + 1;
        else
            i = 1;
        end
        unit = nonClusSet(i);
        putativeClusters = clussMap(find(P(unit,:) == max(P(unit,:))));
        if sum(isnan(putativeClusters)) ~= length(putativeClusters)
           putativeClusters = putativeClusters(~isnan(putativeClusters));
           val = unique(putativeClusters);
           count = histc(putativeClusters, val);
           selCluster = val(find(count == max(count)));
           if length(selCluster) > 1
               dist = [];
               for sC = 1:length(selCluster)
                   dist = [dist norm(sMap.codebook(unit,:) - clusCentroids(sC))];
               end
               selCluster = selCluster(find(dist == min(dist)));
           end
           clussMap(unit) = selCluster;
           nonClusSet(i) = [];
        end
        checkVec = [checkVec length(nonClusSet)];
        %Stop criteria, checks if in the last 100*N iterations the length of
        %nonClusSet has not decreased
        if length(checkVec) > 100*N && length(unique(checkVec(end-(100*N-1):end))) == 1 
            break
        end
    end
    clussMap(find(isnan(clussMap))) = 0;
    clusData = clussMap(BMUs);
    
    % First optional output: number of clusters
    varargout{1} = clusNum;
    
    % Second optional output: association probability matrix
    varargout{2} = W;
    
end