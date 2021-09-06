#' Supervised method to assign cell type
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param method Choosing the published method to assign cell type. Now support method contains: ACTINN, CaSTLe, CEHTAH, Garnett, SciBet, scID, scLearn, scmapcluster, scPred, scVI, Seurat, SingleCellNet, SingleR.
#' @param ... other paramter used in method/classifer, the detial could see the corresponding method.
#'
#' @return The predict cell type.
#'
#' @export
#'
#'
eagersupervised = function(train,
                           test,
                           label_train,
                           method,
                           ...){

  rownames(train) = make.unique(toupper(rownames(train)))
  rownames(test) = make.unique(toupper(rownames(test)))

  common_gene = intersect(rownames(train),rownames(test))

  if(length(common_gene) < 500){
    stop('Please convert the gene name of training dataset and testing dataset into the same format!')
  }

  train = train[common_gene,]
  test = test[common_gene,]

  method = tolower(method)

  predict_label = switch(method,
                         actinn = actinn(train,test,label_train,...),
                         castle = castle(train,test,label_train,...),
                         chetah = chetah(train,test,label_train,...),
                         clustifyr = clustifyr(train,test,label_train,...),
                         garnett = garnett(train,test,label_train,...),
                         schpl = schpl(train,test,label_train,...),
                         scibet = scibet(train,test,label_train,...),
                         scid = scid(train,test,label_train,...),
                         sclearn = sclearn(train,test,label_train,...),
                         scmapcluster = scmapcluster(train,test,label_train,...),
                         scpred = scpred(train,test,label_train,...),
                         scvi = scvi(train,test,label_train,...),
                         seurat = seurat(train,test,label_train,...),
                         singlecellnet = singlecellnet(train,test,label_train,...),
                         singler = singler(train,test,label_train,...),
                         markercount = markercount(train,test,label_train,...),
                         mars = mars(train,test,label_train,...),
                         scclassifr = scclassifr(train,test,label_train,...),
                         stop('Please input valid classifier method, check the method description!'))

  return(predict_label)

}


