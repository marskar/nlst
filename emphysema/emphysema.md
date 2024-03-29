Introduction {#introduction .unnumbered}
============

According to the Global Initiative for Chronic Obstructive Lung Disease
(GOLD), chronic obstructive pulmonary disease (COPD) is a disease
characterized by persistent respiratory symptoms and airflow limitation
caused by a mixture of small airways disease (bronchiolitis) and
parenchymal destruction (emphysema)[@gold19]. COPD is projected to be
the third leading cause of death by 2020 so improved methods for
diagnosis are going to be critical in the coming years. Furthermore, the
primary modality for diagnosis of emphysema is chest CT scan, so there
is substantial potential for quantitative analysis of emphysema.

Quantitative Analysis of Emphysema {#quantitative-analysis-of-emphysema .unnumbered}
----------------------------------

The current state of the art method for performing quantitative analysis
of emphysema is the Density Mask technique[@msma88]. Density Mask
measures the percentage of voxels in the lung below a certain density
threshold, and if this percentage is above a preset threshold then the
patient is considered emphysematous. This technique has shown to be very
effective when used consistenly on the same scanner where the threshold
can be tailored to that scanner and reconstruction kernel; however,
Density Mask suffers greatly when the CT scans are taken from a more
diverse data set as in the case of the NLST[@Gierada10]. The CT scans
from the NLST come from 33 sites across the US and 23 different scanner
models reconstructed using 44 different reconstruction kernels[@nlst11].
Therefore, it is paramount to develop quantitative methods that are much
more robust in the face of ever-growing datasets.

In order to quantitatively analyze emphysema in a more robust manner, we
use as our model a three dimensional convolutional neural network (3D
CNN). Our 3D CNN utilizes multiple layers of convolutions to detect and
weight features at all scales throughout the lung image. This does have
the downside of requiring significantly more computational power, but
most of this is required during training and not during deployment.

Lung Cancer Risk Prediction {#lung-cancer-risk-prediction .unnumbered}
---------------------------

\[\[Cite Hilary’s Markov paper and give a basic explanation here\]\]

Methods {#methods .unnumbered}
=======

The chest CT scans used in this study were taken from the all three
years of the NLST from both the LSS branch and the ACRIN branch; 8416
patients were used for training, 480 for validation, and 2016 for
testing. Each CT scan was labeled by radiologists at the site of the
screening for lung CT abnormalities.

The scans were first converted to NifTI-1 format, then cropped to a
bounding box around the lung using the Progressive Holistically-Nested
Network (P-HNN) lung segmenter[@phnn17], normalized in three different
lung windows of \[-1000, 200\], \[-160, 240\], and \[-1000, -775\], and
rescaled to a standard size of 128x128x128. These lung windows were
chosen due to their use in the P-HNN segmenter[@phnn17]. The resulting
3-channel image was then fed into a standard 3D convolutional neural
network. The network consisted of five 4-layer blocks of 3x3x3
convolution, batch normalization, ReLU activation, and 3D max pooling;
then a convolution group of 2x2x2 convolution, batch normalization, ReLU
activation, and 50% dropout before a fully connected group of 1x1x1
convolution, 50% dropout, 1x1x1 convolution, 50% dropout, a flattening
layer, and a dense layer with 2 class outputs. In order to compensate
for the asymmetry of the labels (there were many more non-emphysematous
cases than emphysematous cases), positive labels were weighted 3 times
as much as negative labels in the training process.

Three neural networks were trained for this experiment. One for just T0
scans, another for T0 and T1 scans, and a final model for all T0, T1,
and T2 scans. Each neural network was trained concurrently on 4 NVIDIA
Titan X graphics cards using Python 2.7 and Keras bindings for
Tensorflow[@keras15]. The majority of the time spent training the model
was spent preprocessing the CT images into the format necessary for our
model.

Results {#results .unnumbered}
=======

As a classifier, the network had moderate success. After 26 epochs of
training, the T0 model had a classification AUC of 0.673; after 20
epochs of training, the T1 model had a classification AUC of 0.689; and
after 24 epochs of training, the T2 model had a classification AUC of
0.684. The confusion matrices for each model are displayed below.

T0 model:

  -- ---------- ---------- ---- ------
                                
     Negative    Positive       
     Negative      1221     39   1260
     Positive      657      99   756
                                
  -- ---------- ---------- ---- ------

T1 model:

  -- ---------- ---------- ----- ------
                                 
     Negative    Positive        
     Negative      1828     579   2407
     Positive      797      844   1641
                                 
  -- ---------- ---------- ----- ------

T2 model:

  -- ---------- ---------- ------ ------
                                  
     Negative    Positive         
     Negative      2797     757    3554
     Positive      1343     1167   2510
                                  
  -- ---------- ---------- ------ ------

Discussion {#discussion .unnumbered}
==========

The results presented here show that 3D CNN-based approaches are a
promising direction of research for computer-aided diagnosis of
emphysema. As well as just diagnosing emphysema, neural networks could
identify multiple CT-diagnosed diseases to identify risk factors for
lung cancer, or even identify lung cancer risk directly.

Conclusions {#conclusions .unnumbered}
===========

While not accurate enough to be used in a clinical environment, this
shows that there is value in working directly with CT data instead of
just in derived values. Furthermore, the neural network approach is
significantly more robust than the density mask approach, even with a
very basic network.
