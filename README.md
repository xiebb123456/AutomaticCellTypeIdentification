# AutomaticCellTypeIdentification

AutomaticCellTypeIdentification is a wrapper of published automatic cell type identification methods which contains supervised methods, unsupervised methods and semi-supervised methods.

<p align="center" width="100%">
    <img src="figure/website.png"> 
</p>


## Installation

You can install AutomaticCellTypeIdentification from github with:

```R
devtools::install_github('xiebb123456/AutomaticCellTypeIdentification')
```
Note: AutomaticCellTypeIdentification is a wrapper of published methods, the needed package is in Description file.

'''docker'''
sudo docker pull registry.cn-hangzhou.aliyuncs.com/xiebb123456/automaticcelltypeidentification
'''

## Running AutomaticCellTypeIdentification methods

Now, three interface of ```supervised```, ```unsupervised```, ```semisupervised``` methods supports the available robotic methods.  
```supervised``` methods include Seurat, ACTINN, CaSTLe, CHETAH, Garnett, SciBet, scID, scLearn, scmapcluster, scPred, scVI, SingleCellNet and SingleR.  
```unsupervised``` methods include CELLBLAST, CellFishing.jl and scmapcell.  
```semisupervised``` methods include SCSA, DigitalCellSorter and SCINA.  

### Prepare input data  
The training data with cell type could download from GEO/ArrayExpress/GSA.  
The canonical marker of cell type could download from PanglaoDB/CellMarker/CancerSEA.  
The training and testing data, the row is gene and the column is cell, the count format is suggested.

### Example with running supervised methods Seurat
```R
supervised(train,test,label_train,method='Seurat')
```

### Example with running unsupervised methods CELLBLAST
```R
unsupervised(train,test,label_train,method='CELLBLAST',python_link='/home/anaconda3/envs/cellblast/bin/python')
```

### Example with running supervised methods SCSA
```R
semisupervised(test,marker,method='SCSA',python_link='/home/anaconda3/envs/scsa/bin/python')
```

## Tutorial
For more details and basic usage see following tutorials:
[Guided Tutorial](vignettes/introduction.Rmd)

## Contact
Feel free to submit an issue or contact us at xiebb7@mail.sysu.edu.cn for problems about the package installation and usage.
