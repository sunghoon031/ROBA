# ROBA: Rotation-Only Bundle Adjustment

[Paper](https://arxiv.org/abs/2011.11724), [Supplementary Material](https://seonghun-lee.github.io/pdf/SupplementaryMaterial_RotationOnlyBundleAdjustment.pdf)

In this repository, we provide the implementation of ROBA. If you use our code, please cite it as follows:

````
@InProceedings{Lee_2021_CVPR,
author = {Lee, Seong Hun and Civera, Javier},
title = {Rotation-Only Bundle Adjustment},
booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {June},
year = {2021}
}
````

### Quick start
1. To initialize the rotations, ROBA relies on rotation averaging. Specifically, we used the rotation averaging code provided by Chatterjee and Govindu. Go to [their webpage](http://www.ee.iisc.ac.in/labs/cvl/research/rotaveraging/) and download their code (PAMI version). Extract `SO3GraphAveraging` folder and save it in the same directory.

2. Download the processed input/results data from [here](https://drive.google.com/drive/folders/1-Bnbk2P8NIfycUJ3_sMY8_D-KugxkSp1?usp=sharing) and save it in the same directory. 

3. In Matlab, run the following command to compile the mex function: `mex computeEdgeCost.cpp`.

4. Run `RunSynthetic.m` (or `RunReal.m`) to try it on the synthetic (or real) data.

### File descriptions

- `RunSynthetic.m`: 
Run the experiment on synthetic data. 
Set `draw3D = true` to visualize the ground truth data.
Choose the simulation configuration by setting the `config` variable.

- `RunReal.m`: 
Run the experiment on the real data.
Set `plotResult = true` to show the accuracy evolution over iterations.
Set `datasets` to select the datasets to evaluate. 
The Booleans `square_rooting`, `switch_alpha` and `approximate_gradient` are for the ablation study in the paper. 

- `computeEdgeCost.cpp`:
C++ Mex implementation of Algorithm 1 in our paper.

- `GetGroundTruthRealData.m`: 
Parse data from [1DSfM datasets](https://www.cs.cornell.edu/projects/1dsfm/). 
It is not necessary to run this script, as we already provide the parsed data [here](https://drive.google.com/drive/folders/1-Bnbk2P8NIfycUJ3_sMY8_D-KugxkSp1?usp=sharing).
If you want to try this yourself, then you need to download the raw data (`.out`, `coords.txt`, `EGs.txt`, `tracks.txt`) from their webpage.

- `LoadAndPlotResultsInPaper.m`: 
Load and plot the published results in our paper. Our results are available [here](https://drive.google.com/drive/folders/1-Bnbk2P8NIfycUJ3_sMY8_D-KugxkSp1?usp=sharing).
