function RhythSOM_plotTaggedEvents(Data, Tags, TagsID, sMap, ClusData, clussMap, BMUs, showNonClus, waveletParams)

%   RhythSOM_plotTaggedEvents plots as many Self-Organising Maps as type of tags.
% Depending on the type of clusterization chosen in "probClus" (in the main script), 
% the figure shows:
%
%    -  If probClus was 1: there is just one figure displaying over each node of 
%       each of these SOM figures is written the number of events 
%       belonging to that node.
%
%    -  If probClus was 0: there are two figures. The first one is the same as described
%       above. The second one plots the clusterized SOM with the mean wave in the middle.
%
% 
%   RhythSOM_plotTaggedEvents(Data, Tags, TagsID, sMap, ClusData, clussMap, BMUs, showNonClus, waveletParams)
% 
% Inputs
%
%
%   Data     (matrix)  N x M matrix, formed by N samples of size M. 
%                      Examples of usable data could be N 100ms-long
%                      somatic sharp wave ripples, or N theta cycles.
%
%   Tags     (matrix)  N x 1 matrix, tags going from 1 up to the number
%                      of tags.
%                      Eg.: Tags' = 
%                               1   3   1   2   1   3   3   1   3
%
%   TagsID (cell array)  1 x #tags cell array, with strings containing the
%                      name of each tag
%                      Eg.: TagsID = 
%                               {'Type1', 'Type2', 'Type3'}
%
%   sMap     (struct)  SOM ouput sctructure from som_make.
%
%   clusData (matrix)  1 x N matrix, N being the number of samples.
%                      For each sample, the cluster number to which
%                      it belongs. 
%                      Eg.: clusData = 
%                               1   1   1   2   2   3   3   3   3
%
%   clusterMap (matrix)  1 x U matrix, U being the number of units on 
%                      the sMap. For each unit, the cluster number to
%                      which it belongs.
%                      Eg.: clussMap = 
%                                   1   1   2   3   3
%
%   BMUs     (matrix)  Best-matching units output matrix from som_bmus
%                      for a given sMap.
%
%   showNonClus (bool) If probClus==1, and if showNonClus==1, the number of
%                      events in black cells are shown, and if is 0, not.
%
%   waveletParams (struct)  Optional. Structure with parameters to plot the 
%                      wavelet on Figure 1: sample frecuency first, 
%                      frequency limits second, and number of consecutive
%                      waveforms to be displayed. If not given,
%                      wavelets are not shown, sample frequency is
%                      taken as 1Hz, and just one repetition of the
%                      event is displayed.
%                      Eg.: 20000 Hz, frequencies between 15 and 250 Hz
%                           and two repeated theta cycles displayed
%
%                              waveletParams = {20000, [15 250], 2}
%
%                      Eg.: 5000 Hz, frequencies between 70 and 500 Hz
%                           and ripple is displayed just once
%  
%                              waveletParams = {5000, [70 500], 1}
%
% 
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%                          
% LCN - E.F Rodríguez Sebastián 2019


    if nargin < 8 || isempty(showNonClus)
        showNonClus = zeros(1,length(TagsID));
    end

    if nargin < 9 || isempty(waveletParams)
        waveletParams = {};
        fs = 1;
        numEventPlot = 1;
    end

    % Parameters for wavelet
    if ~isempty(waveletParams)
        fs = waveletParams{1};
        freqLims = waveletParams{2};
        numEventPlot = waveletParams{3};
    end
    
    
    %Analysis evaluation
    headers = {'ClusterID'};
    for h = 1:length(TagsID)
        headers = [headers strcat(TagsID{h}, 'Number')];
    end
    for h = 1:length(TagsID)
        headers = [headers strcat(TagsID{h}, 'Percentage')];
    end
    percTable = array2table(NaN(max(ClusData),length(headers)),'VariableNames', headers);
    for clus = 1:max(ClusData)
        percTable.ClusterID(clus) = clus;
        data = Data(find(ClusData == clus));
        tags = Tags(find(ClusData == clus));
        total = size(data, 1);
        for t = 1:length(TagsID)
            percTable.([TagsID{t} 'Number'])(clus) = length(find(tags == t));
            percTable.([TagsID{t} 'Percentage'])(clus) = (length(find(tags == t))/total)*100;   
        end
    end
    
    %SOM figure displaying waves and divided in event types (Only for non
    %probabilistic clusters)
    if ~ismember(0, clussMap)
        figure
        for t = 1:length(TagsID)
            data = Data(find(Tags == t),:);
            clusData = ClusData(find(Tags == t));
            if length(TagsID) >= 3 
                subplot(ceil(length(TagsID)/3),3,t)
            else
                subplot(ceil(length(TagsID)/2),2,t)
            end
            title([TagsID{t} ' events'])
            hold on    
            % Cluster colouring
            colors = colormap('hsv');
            Code = colors(round(clussMap/max(clussMap)*length(colors)),:);
            % Plot sMap
            som_cplane(sMap,Code);
            % Plot mean waveform on each cluster or number of events
            clussMapMat = reshape(clussMap,sMap.topol.msize);
            % Width and hight
            W = 0;
            H = 0;
            for iClus = 1:max(ClusData)
                [xs, ys] = find(clussMapMat'==iClus);
                X = mean(xs);
                Y = mean(ys);
                W = W + (max(xs)-min(xs))/max(ClusData);
                H = H + 0.5*(max(ys)-min(ys))/max(ClusData);
            end
            for iClus = 1:max(ClusData)
                [xs, ys] = find(clussMapMat'==iClus);
                X = mean(xs);
                Y = mean(ys);

                xmax = size(Data,2);
                xnorm = [1:xmax]/xmax;
                tempData = data(clusData==iClus,:);
                dataS = std(tempData, 0, 1)/(max(max(data))-min(min(data)));
                dataM = mean(tempData, 1)/(max(max(data))-min(min(data)));
                dataM = dataM - mean(dataM);
                fill( X-W/2+[xnorm,flip(xnorm)]*W, Y+([dataM+dataS,flip(dataM-dataS)])*H,1,'facecolor','w','facealpha',0.6,'edgecolor',[1,1,1]*0.4);
                plot( X-W/2+xnorm*W, Y+dataM*H,'k','linewidth',2);
                text( X-W/4, Y-H*2/4, sprintf('Cluster #%d',iClus) )
                perc = percTable.([TagsID{t} 'Percentage'])(iClus);
                text( X-W/4+0.2, Y-H*2/4+0.3, sprintf('%2.2f %%',perc) )
            end
            ax = gca; outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1)/2;
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height-0.1];
        end
    end
    
    
    %SOM figure displaying number of events and divided in event types
    f = figure('pos',[150 200 700 1000]);
    for t = 1:length(TagsID)
        data = Data(find(Tags == t),:);
        clusData = ClusData(find(Tags == t));
        bmus = BMUs(find(Tags == t));
        if length(TagsID) >= 3 
            subplot(ceil(length(TagsID)/3),3,t)
        else
            subplot(ceil(length(TagsID)/2),2,t)
        end
        title([TagsID{t} ' events (' num2str(size(data,1)) ')'])
        hold on    
        % Cluster colouring
        colors = colormap('hsv');
        if ismember(0, clussMap)
            Code = zeros(length(clussMap), 3);
            values = round(clussMap/max(clussMap)*length(colors));
            for v = 1:length(values)
                if values(v) ~= 0
                    Code(v,:) = colors(values(v),:);
                end
            end
        else
            Code = colors(round(clussMap/max(clussMap)*length(colors)),:);
        end
        % Plot sMap
        som_cplane(sMap,Code);
        % Plot mean waveform on each cluster or number of events
        numEvents = NaN(length(clussMap), 1);
        for u = 1:length(clussMap)
            numEvents(u) = length(find(bmus == u));
        end
        numEventsMat = reshape(numEvents,sMap.topol.msize);
        clussMapMat = reshape(clussMap,sMap.topol.msize);
        for ii = 1:size(numEventsMat, 2)
            for jj = 1:size(numEventsMat, 1)
                if numEventsMat(jj,ii) ~= 0
                    if clussMapMat(jj, ii) ~= 0
                        txt = text(ii, jj, num2str(numEventsMat(jj,ii)), 'FontSize',14-(ceil(length(TagsID)/2)*2 - 2));
                        set(txt, 'HorizontalAlignment', 'center');
                    else
                        if showNonClus(t) == 1
                            txt = text(ii, jj, num2str(numEventsMat(jj,ii)), 'Color', [1 1 1], 'FontSize',14-(ceil(length(TagsID)/2)*2 - 2));
                            set(txt, 'HorizontalAlignment', 'center');
                        end
                    end
                end
            end
        end
    end  
    
   
    %Figures of cluster waves for each type of event
    for t = 1:length(TagsID)
        data = Data(find(Tags == t),:);
        clusData = ClusData(find(Tags == t));
        
        fig = figure('pos',[150 200 700 900]);
        % Subplot numbers in x and y
        spy = round(sqrt(max(ClusData)));
        spx = ceil(max(ClusData)/spy);
        times = [1:size(data,2)]/fs*1000;
        tend = times(end);
        for ii = 2:numEventPlot
            times = [times [1:size(data,2)]/fs*1000+(ii-1)*tend];
        end
        
        for iClus = 1:max(ClusData)
            % Subplot
            subplot(spx,spy,iClus); hold on;
            % Get data to plot
            tempData = data(clusData==iClus,:);
            for ii = 2:numEventPlot
                tempData = [tempData data(clusData==iClus,:)];
            end
            dataS = std(tempData, 0, 1);
            dataM = mean(tempData, 1);

            % Background wavelet
            if ~isempty(waveletParams)
                varNorm = 1./(max(max(data))-min(min(data)))*diff(freqLims)*0.5;
                dt = 1/fs;
                numVoices = 32;
                f0 = centfrq('morl');
                scales = helperCWTTimeFreqVector( freqLims(1), freqLims(2), f0, dt, numVoices);
                cwt = cwtft( {dataM, dt} ,'wavelet','morl','scales',scales);
                helperCWTTimeFreqPlot(cwt.cfs, times, cwt.frequencies, 'surf',[],'Time (ms)','Freq (Hz)');
                ch = colorbar; delete(ch);
                % Mean waveform
                plot3( times, dataM*varNorm + mean(freqLims)*0.9, (max(get(gca,'ZLim'))+0.1)*ones(1,numel(times)), 'k', 'linewidth', 2);
            else
                % Mean +/- std of waveform
                fill( [times,flip(times)], [dataM-dataS, flip(dataM+dataS)], 1,'facecolor',0.7*[1,1,1],'facealpha',0.6,'edgecolor','none')
                % Mean waveform
                plot( times, dataM, 'k', 'linewidth', 2);
            end
            % Axes
            perc = percTable.([TagsID{t} 'Percentage'])(iClus);
            title(sprintf('Cluster #%d (%2.2f %%)',iClus, perc))
        end
        sgtitle([TagsID{t} ' events'])
    end

    
    
    
end
