# Project 2: Machine Learning

In this project, you will write a proposal for your bioinformatics/ML project

> *“Wow, my model has 97% accuracy!” - a happy moment just before realizing the test set had unintentionally leaked into the training data.*

---

## Project Proposal

Your proposal should include the following sections.

### 1. Project Overview

- **Overarching goal**  
  Development of a deep-learning–based tool to resolve batch effects in flow and mass cytometry

- **Rationale**  
  Flow and mass cytometry are widely used single-cell assays, but the measured fluorescence intensity depends on the specific cytometer and fluorochromes used. As a result, ML or DL models trained on data from a single site often fail to generalize to other datasets because of batch effects. If this batch-effect problem can be resolved, it would become possible to develop universal models that apply across cytometers, laboratories, and experimental conditions, thereby streamlining flow cytometry analysis such as autogating.

- **Specific aims** (at least two)  
  1)  
  - **Aim 1**: Develop a model that mitigates batch effects in cytometry data while minimizing loss of underlying biological information.  
  - **Expected outcome**: Substantial reduction of batch effects across FCS files, enabling more consistent cross-batch analysis.  
  - ** Potential Challenge**: Balancing batch-effect correction against preservation of biological signal, and avoiding over-smoothing or loss of rare populations.  
  2)  
  - **Aim 2**:   
  - **Expected outcomes**: Robust integration of flow cytometry data regardless of marker panel differences.  
  - **Potential challenges**: Accurately predicting missing marker expression from a limited set of partially non-orthogonal observed markers.  


### 2. Datac

Describe the dataset you will work with and how you plan to manage it.

- **Dataset description**
  - **Source**: Flow cytometer dataset from Openrepository and immport.
  - **Size**: ~1k of fcs files, 10~20 features, 1 mil ~ 1 bil of cells per each fcs file. 
  - **Format**: .fcs format

- **Data suitability**
  - Should be converted to easy extractable matrix format. Raw fcs file should be reformatted to tabular format like hdf5

- **Storage and data management**
  - Where will you store the dataset?: GCP platform
  - How will you back it up?: rclone
  - How will you share it with collaborators if needed?: Granting access to GCP


### 3. Environment

Document how your computational environment will be set up.

- **Coding environment**
  - Local machine

- **Dependencies**
  - Flowkiy, Pytorch, Numpy and etc.

- **Reproducibility**
  - How will you ensure others can rerun your analysis?
    - Version control with environment.yaml

### 4. Pipeline

Describe the sequence of steps your analysis will follow.

- **Algorithms and methods**  
  What models, algorithms, or computational steps do you plan to run? Are there steps that depends on output of other steps?

- **Scalability and efficiency**  
  How will you ensure your pipeline runs efficiently on your dataset size, format, number of samples?


### 5. Machine Learning

Brainstorm an ML task that can be performed on your data

- **Task definition**  
  For batch-effect mitigation training model should be trained without supervision while for feature imputation model should be trained with supervision. 

- **Feature representation**  
  Tabular dataset with rows representing each cell and column representing each marker. 

- **Model selection**  
  Masked-Graph variable autoencoder. 

- **Generalization strategy** (for supervised learning)  
  Using optimal transport theorem for feature matching 

- **Evaluation metrics**  
  Batch effect evaluation marker like iLISI, Silhouette label, Leiden ARI, Leiden NMI, KBET, Graph connectivity and etc... 
