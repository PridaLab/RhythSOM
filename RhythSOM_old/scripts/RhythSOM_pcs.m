function DataPCA = RhythSOM_pcs(Data, minPCvar)

%   RhythSOM_pcs computes a Principal Component Analysis (PCA), and 
% returns the first Principal Components (PCs). The number of PCs can
% be selected manually (if minPCvar is 'nan'), or automatically (if 
% minPCvar is a number between 0 and 1).
% 
% DataPCA = RhythSOM_pcs(Data, minPCvar)
% 
% Inputs
%
%   Data     (matrix)  N x M matrix, formed by N samples of size M. 
%					   Examples of usable data could be N 100ms-long
% 					   somatic sharp wave ripples, or N theta cycles.
%
%   minPCvar (num/nan) Optional. Minimum explained variance for the
%                      selected number of Principal Components. For
%                      example, if 'minPCvar=0.99',the number of PCs
%                      will be chosen automatically to be the one that
%                      explains 99% of variance. If 'minPCvar=nan',
%                      then the number of PCs will be chosen manually
%                      based on a graph. Variable can be then 'nan' or
%                      a float between 0 and 1.
%
% 
% Outputs
%
%   DataPCA  (matrix)  N x K matrix, formed by N samples of size K 
%                      (selected number of PCs). It's the Principal
%                      Components matrix of Data.
%                      
% LCN-acnavasolive 2019


    % Principal component analysis
    [~ , DataPCA , sumVariance] = pca(Data);

    % Select number of principal components as closest to minPCvar
    [~, numPCs] = min( abs( minPCvar - cumsum(sumVariance)/max(cumsum(sumVariance)) ) );
    
    % Automatic selection
    if ~isnan(minPCvar)
        % Figure
        fig = figure('position',[100,300,500,400]);
        grid on; hold on;
        % Explained variance for different total number of PCs
        semilogx(cumsum(sumVariance)/max(cumsum(sumVariance)),'k','linewidth',2);
        plot([ numPCs numPCs ], [0 1.01],'r','linewidth',2);
        xlabel('PCs'); ylabel('PC variance (cummulative)')
        title('Selected number of PCs');
        xtickVec = [1, 10, 25, 50, 100, 250, 500, 1000, 2000 5000];
        xtickVec = xtickVec( xtickVec < length(sumVariance) );
        set(gca,'xtick',xtickVec)
        axis([0 length(sumVariance) 0 1.01])
        fprintf('   - Number of PCs: %d\n',numPCs);
    
    % Manual selection
    else
        % Selecting bar
        hbar = [numPCs numPCs];
        % Figure
        figure('position',[100,300,500,400]); grid on; hold on;
        % Variation of each PC
        semilogx(cumsum(sumVariance)/max(cumsum(sumVariance)),'k','linewidth',2);
        xlabel('PCs'); ylabel('PC variance (cummulative)')
        title('Select number of PCs (drag bar)');
        xtickVec = [1, 10, 25, 50, 100, 250, 500, 1000, 2000 5000];
        xtickVec = xtickVec( xtickVec < length(sumVariance) );
        set(gca,'xtick',xtickVec)
        axis([0 length(sumVariance) 0 1.01])
        % Bar to select number of PCs
        ax = axis;
        plotPCs = plot(hbar, [ax(1) ax(2)],'r','linewidth',2);
        draggable(plotPCs,'constraint','h');
        get(plotPCs,'XData');
        btn = uicontrol('Style', 'pushbutton', 'String', 'Select',...
                        'Units','normalize','Position', [.8 .94 .15 .06],...
                        'Callback', 'uiresume(gcbf)');
        uiwait(gcf);
        numPCs = round(mean(get(plotPCs,'XData')));
        fprintf('   - Number of PCs: %d\n',numPCs);
    end
    
    % Then take the first <numPCs> number of PCs
    DataPCA = DataPCA( : , 1:numPCs );
    
    
end