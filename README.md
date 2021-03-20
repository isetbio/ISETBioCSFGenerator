## ISETBio CSF Generator
Flexiable framework for simulation of contrast sensitivity function and threshold related visual task in ISETBio.

## Getting Started
We have written a set of well-commented tutorial code (under `tutorials`) to demonstrate the basic usage of the CSF generator.  

- `t_spatialCSF.m` uses the CSF Generator to run out CSFs in different color directions.  
- `t_chromaticThresholdContour.m` computes isothreshold contour in different color directions.

## Design
The CSF generator is consist of four main components:
- Scene Engine: Procedure for stimulus generation (see `t_sceneGeneration.m`). 
- Neural Response Engine: Produces neural responses given a visual stimulus (see `t_neuralResponseCompute.m`).  
- Response Classifier Engine: Produces binary decisions based on neural responses (see `t_responseClassifier.m`).
- Quest Threshold Engine: Adaptive procedure for measuring psychometric curve (see `t_thresholdEngine.m`).

Please check `toolbox` for the interface convention for each of the class.  

