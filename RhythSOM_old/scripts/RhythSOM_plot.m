function RhythSOM_plot(Data, sMap, clusData, clussMap, waveletParams, dirSave)

%   RhythSOM_plot plots to images that show the waveform clusterization.
%   First image (on the left) shows the mean waveform in black along with
% its standard deviation in gray. If the waveletParams input is given, 
% the wavelet of this mean waveform is plotted on the background.
%   Second image (on the right) shows the unfolded self-organized map.
% Each colored area belongs to a different cluster, and areas near on the
% map indicate similarity. Over each coloured area, the mean waveform is
% again plotted.
% 
%   RhythSOM_plot(Data, sMap, clusData, clussMap, waveletParams, dirSave)
% 
% Inputs
%
%   Data         (matrix)  N x M matrix, formed by N samples of size M. 
%					       Examples of usable data could be N 100ms-long
% 					       somatic sharp wave ripples, or N theta cycles.
%
%   sMap          (struct)  SOM ouput sctructure from som_make.
%
%   clusData      (matrix)  1 x N matrix, N being the number of samples.
%                           For each sample, the cluster number to which
%                           it belongs. 
%                           Eg.: clusData = 
%                                    1   1   1   2   2   3   3   3   3
%
%   clussMap      (matrix)  1 x U matrix, U being the number of units on 
%                           the sMap. For each unit, the cluster number to
%                           which it belongs.
%                           Eg.: clussMap = 
%                                    1   1   2   3   3
%
%   waveletParams (struct)  Optional. Structure with parameters to plot the 
%                           wavelet on Figure 1: sample frecuency first, 
%                           frequency limits second, and number of consecutive
%                           waveforms to be displayed. If not given,
%                           wavelets are not shown, sample frequency is
%                           taken as 1Hz, and just one repetition of the
%                           event is displayed.
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
%   Figure 1     (figure)  'RhythSOM_meanWaveforms.png'. Save on dirSave if
%                          given; if not, in the current folder.
%
%   Figure 2     (figure)  'RhythSOM_map.png'. Save on dirSave if given;
%                          if not, in the current folder.
%
%
% Requires SOM Toolbox 2.0. Available at: 
%
%    http://www.cis.hut.fi/projects/somtoolbox/
%                          
% LCN-acnavasolive 2019

    if nargin < 4; dirSave = ''; end
    if nargin < 3 || isempty(waveletParams)
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
    % Number of clusters
    clusNum = max(clusData);
    
    
    % FIGURE 1
    fig1 = figure('pos',[150 200 700 900]);
    % Subplot numbers in x and y
    spy = round(sqrt(clusNum));
    spx = ceil(clusNum/spy);
    times = [1:size(Data,2)]/fs*1000;
    tend = times(end);
    for ii = 2:numEventPlot
        times = [times [1:size(Data,2)]/fs*1000+(ii-1)*tend];
    end
    for iClus = 1:clusNum
        % Subplot
        subplot(spx,spy,iClus); hold on;
        % Get data to plot
        data = Data(clusData==iClus,:);
        for ii = 2:numEventPlot
            data = [data Data(clusData==iClus,:)];
        end
        dataS = std(data);
        dataM = mean(data);
        
        % Background wavelet
        if ~isempty(waveletParams)
            varNorm = 1./(max(max(Data))-min(min(Data)))*diff(freqLims)*0.5;
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
        title(sprintf('Cluster #%d',iClus))
    end
    
    
    % FIGURE 2    
    fig2 = figure('pos',[950 200 700 900]); hold on
    % Cluster colouring
    colors = colormap('hsv');
    Code = colors(round(clussMap/max(clussMap)*length(colors)),:);
    % Plot sMap
    som_cplane(sMap,Code);
    % Plot mean waveform on each cluster
    clussMapMat = reshape(clussMap,sMap.topol.msize);
    % Width and hight
    W = 0;
    H = 0;
    for iClus = 1:clusNum
        [xs, ys] = find(clussMapMat'==iClus);
        X = mean(xs);
        Y = mean(ys);
        W = W + (max(xs)-min(xs))/clusNum;
        H = H + 0.5*(max(ys)-min(ys))/clusNum;
    end
    for iClus = 1:clusNum
        [xs, ys] = find(clussMapMat'==iClus);
        X = mean(xs);
        Y = mean(ys);
        
        xmax = size(Data,2);
        xnorm = [1:xmax]/xmax;
        data = Data(clusData==iClus,:);
        dataS = std(data)/(max(max(Data))-min(min(Data)));
        dataM = mean(data)/(max(max(Data))-min(min(Data)));
        dataM = dataM - mean(dataM);
        fill( X-W/2+[xnorm,flip(xnorm)]*W, Y+([dataM+dataS,flip(dataM-dataS)])*H,1,'facecolor','w','facealpha',0.6,'edgecolor',[1,1,1]*0.4);
        plot( X-W/2+xnorm*W, Y+dataM*H,'k','linewidth',2);
        text( X-W/8, Y-H*3/4, sprintf('Cluster #%d',iClus) )
    end
    ax = gca; outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1)/2;
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height]; 
    
    % Save images
    saveas(fig1,[dirSave,'RhythSOM_meanWaveforms.png']);
    saveas(fig2,[dirSave,'RhythSOM_map.png'])

end















