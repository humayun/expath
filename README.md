# ExPath
ExPath Image Classification Framework

Expansion Pathology Iamge Classification: This framework classify Expansion pathology images into Normal, Benign Breast Lesions (UDH and ADH) and Pre-invasvie Breast Leasions. This framework consisted of three components:

1- Nuclei Segmentation and Feature Extraction written in matlab

2- Nuclei Segmentation Evalaution written in matlab using LabelMe API

3- Image Classification written in R

This image classification framework for post-expansion DAPI-stained images includes foreground detection, nucleus seed detection, and nuclear segmentation. Following application of the framework, we extract three kinds of features from each segmented nucleus from both the pre-expanded and post-expanded images: nuclear morphology features, nuclear intensity features, and nuclear texture features.

