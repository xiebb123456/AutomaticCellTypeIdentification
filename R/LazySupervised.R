#' UnSupervised method to assign cell type
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param method Choosing the published method to assign cell type. Now support method contains: CELLBLAST, CellFishing.ji,scmapcell
#' @param ... other paramter used in method/classifer, the detial could see the corresponding method.
#'
#' @return The predict cell type.
#'
#' @export
#'
#'
lazysupervised = function(train,
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
                         cellblast = cellblast(train,test,label_train,...),
                         cellfishing.ji = cellfishing.ji(train,test,label_train,...),
                         scmapcell = scmapcell(train,test,label_train,...),
                         stop('Please input valid classifier method, check the method description!'))

  return(predict_label)

}






