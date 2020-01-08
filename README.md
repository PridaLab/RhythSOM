# RhythSOM

RhythSOM is a computational tool for MATLAB thought to be used in electrophysiological labs to classify brain rhythms or events, such hippocampal sharp wave ripples or theta cycles, in a unsupervised way, just by their waveform.

## Prerequisites

RhythSOM is based on the [SOM Toolbox 2.0](http://www.cis.hut.fi/projects/somtoolbox/)  MATLAB package, so before using it, it is necessary to download it [(click here)](http://www.cis.hut.fi/projects/somtoolbox/download/).

Once downloaded it, and set into MATLAB path (e.g. by downloading it into your project folder, or by using *addpath()* and *genpath()*), you can start using RhythSOM.

Disclaimer: If some modifications of the toolbox are necessary they will be saved in the Toolbox modifications folder.

## Usage

RhythSOM main function is `RhythSOM_Classifier()`. It classifies the input data into different clusters, and returns an array that tags each sample with the number of the cluster it has been classified to. More information about the process is in next section, `The Code > RhythSOM_Classifier.m`.

### Inputs

* `Data`: N x M matrix, formed by N samples of size M, so each row is a sample.
For instance, data input could be `N=5.000` 50ms-long sharp wave ripples from CA1 Stratum Pyramidale, or `N=10.000` 150ms-long theta cycles from CA1 *Stratum Pyramidale* LFP or LFP from CA1 *Stratum Lacunosum Moleculare*.
* `minPCvar`: `float` between 0 and 1, or a `nan` (Optional). It is minimum explained variance for the selected number of Principal Components. For example, if `minPCvar=0.99`,the number of PCs will be chosen automatically to be the one that explains 99% of variance. If `minPCvar=nan`, then the number of PCs will be chosen manually based on a graph.
* `probClus`: boolean (Optional). Executes custom probabilistic version of k-means (1) or classical version of k-means (0). classical version by default.
* `autoClus`: boolean (Optional). Select number of clusters automatically (1) or manually (0)? Automatically by default.
* `numRep`: `int` (Optional). Number of repetitions for k-means algorithm. If numRep is greater than 1, several iterations will be made, and results from best iteration will be selected. This way the performance doesn't depend so much on local minima near the first random centroids.
* `waveletParams`: MATLAB `structure` (Optional). Structure with parameters to plot the wavelet in the output figures. Parameters are: 1) sample frecuency, 2) wavelet frequency limits, and 3) number of consecutive waveforms to be displayed. If `waveletParams` is not given, wavelets are not shown.
    Eg.: 20000 Hz, frequencies between 15 and 250 Hz and two repeated theta cycles displayed
    ```
    waveletParams = {20000, [15 250], 2}
    ```
    Eg.: 5000 Hz, frequencies between 70 and 500 Hz and ripple is displayed just once .
    ```
    waveletParams = {5000, [70 500], 1}
    ```
* `dirSave`: `string` (Optional). Directory to save the ouput figures. If not given, figures are saved on the current location.
    Eg.: '/home/Projects/RhythSOM/images/'

### Outputs

* `clusData`: 1 x N matrix, N being the number of samples. For each sample, the cluster number to which it belongs.
	```
    Eg.: clusData =
            1   1   1   2   2   3   3   3   3
    ```
* clussMap: matrix. 1 x U matrix (Optional). For each unit of the Self-Organized Map, SOM, (U being the number of units), the cluster number to which it belongs.
	```
    Eg.: clussMap =
            1   1   2   3   3
    ```
* `sMap`: `struct` (Optional). SOM ouput sctructure from `som_make`.
* Figure 1: *figure*. Saved as 'RhythSOM_meanWaveforms.png' on `dirSave`; if `dirSave` input not given, in the current folder.
* Figure 2: *figure*. Saved as 'RhythSOM_map.png' on `dirSave`; if `dirSave` input not given, in the current folder.

![Example of *RhythSOM_meanWaveforms* figure](https://github.com/acnavasolive/RhythSOM/images/RhythSOM_meanWaveforms.png)

![Example of *RhythSOM_map* figure](https://github.com/acnavasolive/RhythSOM/images/RhythSOM_map.png)

## Examples of usage

**Example 1**. Only input `Data`, and only output two figures.
```
Data = [ sample1; ...
		 sample2; ...
		 ...
		 sampleN ];
RhythSOM_Classifier(Data);
```
Figures saved in current path, and by default parameters used will be:
* `minPCvar` = 0.8 (PCs will explain 80% of variance)
* `autoClus` = 1 (number of clusters selected automatically)
* `numRep` = 100 (100 k-means realisations)
* `waveletParams` = {} (wavelet unplotted)

**Example 2**. All inputs, just `clusData` output.
```
Data = [ sample1; ...
		 sample2; ...
		 ...
		 sampleN ];
minPCvar = 0.9;
probClus = 0;
autoClus = 0;
numRep = 50;
waveletParams = {20000, [20 150], 2};
dirSave = '/home/Projects/RhythSOM/images/';

clusData = RhythSOM_Classifier(Data, minPCvar, probClus, autoClus, numRep, waveletParams, dirSave);

>>> clusData'

ans =

	1	1	2	2	1	3	... clusNum
```
You will be asked to select the total number of clusters `clusNum` draggin a bar in a figure.



## The Code

There are five principal functions in this module:

### RhythSOM_Classifier.m

**RhythSOM_Classifier** is RhythSOM main function. Computational steps are the following:
* A 'Data' matrix with eachsample in a row first goes through a Principal Components Analysis (PCA)
* The selected PCs are fed into a Self-Organizing Map (SOM) algorithm from the [SOM Toolbox](http://www.cis.hut.fi/projects/somtoolbox/).
* Units from the map are clusterized by a k-means algorithm
* Results of these clusterization are finally plotted in two separate figures. On the left,for each cluster mean waveforms and their respective wavelet, and on theright, the Self-Organized Map with the mean wavelet on top.
  (Depending on the type of k-means executed the final plot displays the same information in a slightly different way)

There are three parameters that can set different options:

* The number of Principal Components can be chosen manually (*minPCvar=nan*), or automatically by selecting the minimum explained variance with *minPCvar*

* The type of k-means executed, classical (*probClus=0*), or custom probabilistic (*probClus=1*).

* The number of clusters can be chosen also manually (*autoClus=0*), or   automatically (*autoClus=1*).

The way to use this is:
```
clusData, [clussMap, sMap] = RhythmSOM_Classifier(Data, [minPCvar, probClus, autoClus, numRep, waveletParams, imageDir])
```
where variables inside brakets [] are optional inputs or outputs.

### RhythSOM_pcs.m

RhythSOM_pcs computes a Principal Component Analysis (PCA), and returns the first Principal Components (PCs). The number of PCs can be selected manually (if *minPCvar* is 'nan'), or automatically (if *minPCvar* is a number between 0 and 1).

```
DataPCA = RhythSOM_pcs(Data, minPCvar)
```

### RhythSOM_normdata.m

RhythSOM_normdata computes normalized DataPCA so SOM algorithm works with normalized euclidean distances, and avoid biases.

```
[DataNorm, DataDenorm] = RhythSOM_normdata(DataPCA)
```

### RhythSOM_clusters.m

RhythSOM_clusters clusterizes the sMap Best Matching Units (BMUs) through the k-means method. A metric of how well clustering has been  done is computed, the Davies-Boulding index (the smaller the better). Through 'autoClus' it can be set to choose manually or automatically the number of clusters. Number of clusters with minimum DB index will chosen in the automatic performance. A DB index graph will be shown otherwise.

```
[clusData, clussMap, clusNum, clusCentroids, clusDBidx] = RhythSOM_clusters(sMap, BMUs, method, autoClus, numRep)
```

### RhythSOM_probabilisticClusters.m

RhythSOM_probabilisticClusters clusterizes the sMap Best Matching Units (BMUs) through a variation of the k-means method. The probability of being in the same cluster is computed for each pair of Units and then they are clusterized using the maximum probabilities.

```
[clusData, clussMap, clusNum, probMatrix] = RhythSOM_clusters(sMap, BMUs, method, numTrials)
```

### RhythSOM_plot.m

RhythSOM_plot plots to images that show the waveform clusterization.
First image (on the left) shows the mean waveform in black along with its standard deviation in gray. If the waveletParams input is given, the wavelet of this mean waveform is plotted on the background.
Second image (on the right) shows the unfolded self-organized map. Each colored area belongs to a different cluster, and areas near on the map indicate similarity. Over each colored area, the mean waveform is again plotted.

![alt text](https://github.com/acnavasolive/RhythSOM/blob/master/images/RhythSOM_plot_output.png)

```
RhythSOM_plot(Data, sMap, clusData, clussMap, waveletParams, dirSave)
```

### RhythSOM_plotTaggedEvents.m

RhythSOM_plotTaggedEvents plots as many Self-Organising Maps as type of tags.
Depending on the type of clusterization chosen in "probClus" (in the main script), 
the figure shows:
   -  If probClus was 1: there is just one figure displaying over each node of 
      each of these SOM figures is written the number of events 
      belonging to that node.
   -  If probClus was 0: there are two figures. The first one is the same as described
      above. The second one plots the clusterized SOM with the mean wave in the middle.

```
RhythSOM_plotTaggedEvents(Data, Tags, TagsID, sMap, ClusData, clussMap, BMUs, showNonClus, waveletParams)
```