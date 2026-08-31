# Single-molecule-Confocal-Images-Processing
This repository implements a MATLAB-based quantitative imaging pipeline for extracting single-molecule mRNA and translational signals from 3-D embryonic microscopy data, assigning these molecular objects to individual nuclei, and mapping their quantitative properties onto embryo-scale spatial coordinates.


## 1. Workflow
**The pipeline is therefore designed to connect subcellular single-molecule measurements with tissue-scale spatial patterning.**

* 3-D image import: DAPI, mRNA, and GFP (translational signals) z-stacks for multichannel spatial analysis.
* Nuclear segmentation
* Transcription-focus exclusion: Transcription foci and nucleus-attached mRNAs are exclude from cytoplasmic single-molecule mRNA detection.
* Single-molecule mRNA detection: A preliminary 2-D threshold determined by identify the elbow-point of a curve showing mask count changes following series of tested threshold values. Multi-z consistency filtering, local maxima identification, and morphological refinement.
* mRNA Gaussian fitting: Fit individual RNA candidates with a 2-D Gaussian PSF model to calculate integrated intensity.
* mRNA quality control and copy-number calibration: Estimate single-mRNA intensity from a constrained multi-Gaussian intensity model, converting mRNA intensity to absolute counts. 
* GFP-foci detection: Similar procedure with mRNA identification, except the absolute count conversion.
* Physical-coordinate conversion: Convert x, y, and z coordinates into nanometres using the anisotropic imaging scale before spatial-distance calculations.
* Molecule-to-nucleus assignment: No cell membranes in *Drosophila* embryo at this developmental stage, implement a Voronoi-like segmentation, assinging each RNA and GFP object to its nearest nuclear centroid using 3-D nearest-neighbor analysis.
* RNA–GFP colocalization: 300-nm distance threshold.
* Single-molecule translation quantification: Combine GFP intensity with estimated RNA copy number to derive translation-associated measurements for individual molecular objects.
* Cross-magnification image registration: Register high-magnification molecular data to whole-embryo-scale reference images and transform molecular and nuclear coordinates accordingly.
* Anterior–posterior coordinate mapping: Extract embryo geometry and convert registered positions into normalized AP coordinates.
* Spatial-profile generation: from single-molecule signals to whole-embryo-scale spatial profiles.


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
