#' Semi-supervised method to assign cell type
#' @param test A gene-cell matrix of single cell expression data.
#' @param marker_list A list contains marker gene of cell type.
#' @param method Choosing the published method to assign cell type.
#' @param ... other paramter used in method/classifer, the detial could see the corresponding method.
#'
#' @return The predict cell type.
#'
#' @export
#'
#'
markersupervised = function(test,
                            marker_list,
                            method,
                            ...){

  method = tolower(method)

  predict_label = switch(method,
                         digitalcellsorter = digitalcellsorter(test,marker_list,...),
                         scina = scina(test,marker_list,...),
                         # cellassign = cellassign_r(test,marker_list),
                         scsa = scsa(test,marker_list,...),
                         sctyper = sctyper(test,marker_list,...),
                         markercount = markercount_semi(test,marker_list,...),
                         stop('Please input valid method, you can see description in detail!'))

  return(predict_label)

}
