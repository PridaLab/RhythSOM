function [clusData, clussMap, varargout] = RhythSOM_clusters(sMap, BMUs, autoClus, numRep)


%   RhythSOM_clusters clusterizes the sMap Best Matching Units (BMUs) 
% through the k-means method. A metric of how well clustering has been 
% done is computed, the Davies-Boulding index (the smaller the better).
% Through 'autoClus' it can be set to choose manually or automatically
% the number of clusters. Number of clusters with minimum DB index will
% chosen in the automatic performance. A DB index graph will be shown
% otherwise.
% 
%   [clusData, clussMap, clusNum, clusCentroids, clusDBidx] = RhythSOM_clusters(sMap, BMUs, autoClus, numRep)
% 
% Inputs
%
%   sMap     (struct)  SOM ouput sctructure from som_make.
%
%   BMUs     (matrix)  Best-matching units output matrix from som_bmus
%                      for a given sMap.
%
%   autoClus   (bool)  Optional. Select number of clusters automatically
%                      (1) or manually (0)? Automatically by default.
%
%   numRep      (int)  Optional. Number of repetitions for k-means
%                      algorithm. If numRep is greater than 1, several
%                      iterations will be made, and results from best
%                      iteration will be selected. This way the performance
%                      doesn't depend so much on local minima near the
%                      first random centroids.
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
%   clussMap     (matrix)  1 x U matrix, U being the number of units on 
%                          the sMap. For each unit, the cluster number to
%                          which it belongs.
%                          Eg.: clussMap = 
%                                   1   1   2   3   3
%
%   clusNum         (int)  Optional. Number of clusters.
%
%   clusCentroids (cells)  Optional. Centroids coordinates for final set.
%
%   clusDBidx     (cells)  Optional. Davies-Bouldin index for final set.
%
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%                          
% LCN-acnavasolive 2019
        
    % Clustering with k-means algorithm
    % - Initialization of variables
    tmpclusCentroids = kmeans_clusters(sMap);
    clusCentroids = cell(size(tmpclusCentroids,1),numRep);
    clusIndex = cell(size(tmpclusCentroids,1),numRep);
    clusError = ones(size(tmpclusCentroids,1),numRep);
    clusDBidx = ones(size(tmpclusCentroids,1),numRep);
    % - Several repetitions
    parfor iRep = 1:numRep
        [tmpclusCentroids, tmpclusIndex, tmpclusError, tmpclusDBidx] = kmeans_clusters(sMap);
        clusCentroids(:,iRep) = tmpclusCentroids;
        clusIndex(:,iRep) = tmpclusIndex;
        clusError(:,iRep) = tmpclusError';
        clusDBidx(:,iRep) = tmpclusDBidx';
    end

    % Best Davies-Boulding index (metric of how well clustering has been done)
    [clusNum, clusNumRep] = find(clusDBidx==min(clusDBidx,[],'all'));
    if length(unique(clusNum))==1
        clusNum = clusNum(1);
        clusNumRep = datasample(clusNumRep,1);
    end
    
    % Plot Davies-Bouldin indexes and select number of clusters
    if autoClus 
        % Figure
        figure('pos',[200,400,700,400]); grid on; hold on;
        % DB indexes
        plot(clusDBidx,'o','markerfacecolor',[1,1,1]*0.6, 'markeredgecolor',[1,1,1]*0.8);
        plot(1:size(clusDBidx,1),mean(clusDBidx,2),'.-k','linewidth',1.5)
        plot(clusNum,clusDBidx(clusNum,clusNumRep),'xr','linewidth',1.5)
        quiver(clusNum,max(clusDBidx(clusNum,:))+max(max(clusDBidx))*0.6,0,-max(max(clusDBidx))*0.5,0,'r','linewidth',2)
        xlabel('Number of clusters');
        ylabel('Davies-Bouldin index')
        title('Selected number of clusters (automatic)');
        axis([1.5 size(clusDBidx,1)+0.5 min(min(clusDBidx))*0.9 max(max(clusDBidx))*1.1])
        fprintf('   - Number of clusters: %d\n',clusNum);        
    else
        % Figure
        figure('pos',[200,400,700,400]); grid on; hold on;
        % DB indexes
        plot(clusDBidx,'o','markerfacecolor',[1,1,1]*0.6, 'markeredgecolor',[1,1,1]*0.8);
        plot(1:size(clusDBidx,1),mean(clusDBidx,2),'.-k','linewidth',1.5)
        plot(clusNum,clusDBidx(clusNum,clusNumRep),'xr','linewidth',1.5)
        xlabel('Number of clusters');
        ylabel('Davies-Bouldin index')
        title('Select number of clusters (drag bar)');
        axis([1.5 size(clusDBidx,1)+0.5 min(min(clusDBidx))*0.9 max(max(clusDBidx))*1.1])
        % Bar to select number of PCs
        ax = axis;
        hbar = [clusNum clusNum];
        plotPCs = plot(hbar, [ax(3) ax(4)],'r','linewidth',2);
        draggable(plotPCs,'constraint','h');
        get(plotPCs,'XData');
        btn = uicontrol('Style', 'pushbutton', 'String', 'Select',...
                        'Units','normalize','Position', [.8 .94 .15 .06],...
                        'Callback', 'uiresume(gcbf)');
        uiwait(gcf);
        clusNum = round(mean(get(plotPCs,'XData')));
        delete(btn);
        [~, clusNumRep] = min(clusDBidx(clusNum,:));
        fprintf('   - Number of clusters: %d\n',clusNum);
    end
    
    
    % Clustering based on 'numClus'
    % - Cluster of each BMU
    clussMap = clusIndex{clusNum,clusNumRep};
    % - Cluster of each Data point 
    clusData = clussMap(BMUs);
    
    % First optional output: number of clusters
    varargout{1} = clusNum;
    % Second optional output: centroid of clusters
    varargout{2} = clusCentroids;
    % Third optional output: Davies-Bouldin indexes
    varargout{3} = clusDBidx;

end