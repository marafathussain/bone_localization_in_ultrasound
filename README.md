# Bone Localization in Ultrasound
Fluoroscopy remains the primary intraoperative imaging modality for bone boundary visualization in computer assisted orthopaedic surgery systems. The associated radiation exposure posing risks to both patients and surgical teams gave rise to recent interest in safer non-ionizing real-time intraoperative imaging alternatives such as US. In US guided surgical intervention, bone localization in US images is essential for visualization and guidance, e.g., during fragment positioning in fracture reduction surgeries. Despite recent advancement in US intensity-based automatic bone segmentation, results remain unpredictable due to the high levels of speckle noise, reverberation, and signal drop out. To tackle this challenge, we proposed an effective way to extract 3D bone surfaces using a surface growing approach that is seeded from 2D bone contours. The initial 2D bone contours are estimated from a combination of ultrasound strain images and envelope power images. This work is divided into 2 parts: (1) 2D bone delineation and (2) 3D bone surface generation. 

## Dependency: 
Please download these codes [AM2D](https://drive.google.com/drive/folders/1YvwJpDvTbegXX38hevTk-bX1b79h9k7C?usp=sharing) (Thanks to Prof. Hassan Rivaz), [ZPclustering](https://drive.google.com/drive/folders/19VYSLNLCDRiZhjlSNBvysQW-GDL2H-31?usp=sharing), and add their paths to the code.

## (1) Delineation of 2D bone surface in ultrasound using combined elastography and envelope signal power
We estimate the 2D bone contour from a pair of RF frames taken with and without a little compression. To estimate the bone surface location in this pair of 2D images, we used a real-time strain imaging technique developed by Rivaz et al. (2011), which is based on an analytic minimization of a regularized cost function. The resulting strain image was fused with the envelope power map using a weight. The weight was selected based on an empirical analysis of the mean absolute error (MAE) between the actual and estimated bone surfaces for different weight values. The weight for which the MAE was lowest in this pilot data set was used throughout our experiments. Then we used local linear fits over the maximum intensity point along each scan line of the fused map to produce the final bone boundary. Our initial work has been published in this MICCAI 2014 [paper](https://arafathm.github.io/assets/pdf/MICCAI2014.pdf). To use this method, run [2D_miccai14.m](https://github.com/marafathussain/bone_localization_in_ultrasound/blob/main/2D_miccai14.m). We later published an improved version of this work in Ultrasound Medicine and Biology [paper](https://arafathm.github.io/assets/pdf/Arafat_UMB17.pdf). To use this method, run [2D_umb16.m](https://github.com/marafathussain/bone_localization_in_ultrasound/blob/main/2D_umb16.m).   

## (2) Detection of 3D bone surface in ultrasound seeded from 2D strain
After having identified a seed contour in the 2D elastographic image (1. above), we next use a surface growing approach to extend this contour laterally through the 3D volume set. We published this work in Ultrasound Medicine and Biology [paper](https://arafathm.github.io/assets/pdf/Arafat_UMB17.pdf). To use this method, run [3D_umb16.m](https://github.com/marafathussain/bone_localization_in_ultrasound/blob/main/3D_umb16.m).   

## Citations
If you find these works useful, please cite following papers:
```
@inproceedings{hussain2014robust,
  title={Robust bone detection in ultrasound using combined strain imaging and envelope signal power detection},
  author={Hussain, Mohammad Arafat and Hodgson, Antony and Abugharbieh, Rafeef},
  booktitle={International conference on medical image computing and computer-assisted intervention},
  pages={356--363},
  year={2014},
  organization={Springer}
}
```

```
@article{hussain2017strain,
  title={Strain-initialized robust bone surface detection in 3-D ultrasound},
  author={Hussain, Mohammad Arafat and Hodgson, Antony J and Abugharbieh, Rafeef},
  journal={Ultrasound in Medicine \& Biology},
  volume={43},
  number={3},
  pages={648--661},
  year={2017},
  publisher={Elsevier}
}
```
}
