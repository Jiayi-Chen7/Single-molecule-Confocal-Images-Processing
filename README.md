# Single-molecule-Confocal-Images-Processing
This repository contains codes for image processing of single-molecule mRNAs and corresponding translational signals in 3D embryo volume across spatial resolution scales, related to manuscript "mRNA concentration–dependent translation enables rapid and sharp patterning in resource-constrained Drosophila embryos".

## 1. Basic procedure
* 3D segmentation of both nuclei, single-molecule mRNA, and translational signals along anterior-posterior axis of the Drosophila embryo.

## 2. Running the codes

* Add the related "Function" folder to the search path of current session.
* Run `.m` files with "Main" in filename.
* Add input images in the same path of "Main" file or add the folder containing input images to the search path.
* Codes are MATLAB files (MATLAB R2025a) and can be run on MATLAB software.

## 3. Citation
Please cite the following paper:
Chen, Jiayi, et al. "mRNA concentration–dependent translation enables rapid and sharp patterning in resource constraint Drosophila embryos." bioRxiv (2026): 2026-06.

## 4. Reference
*MultiLayerSpotIdentify.m is modified from spmask.m from 
Xu, Heng, et al. "Combining protein and mRNA quantification to decipher transcriptional regulation." Nature methods 12.8 (2015): 739-742.

*GetElbowIndex.m is modified from Jona's answer on Stackoverflow
https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

*parameter filters after PSF fitting are adapted from spfilter.m from
Xu, Heng, et al. "Combining protein and mRNA quantification to decipher transcriptional regulation." Nature methods 12.8 (2015): 739-742.

## 5. Contact Information
- Jiayi Chen (jiayi.chen@pku.edu.cn)
- Feng Liu (liufeng@hebut.edu.cn)
