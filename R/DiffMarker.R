#' Differential marker gene of cell type
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param label_train A vector of cell type label of train.
#' @param only.pos Only return positive markers (FALSE by default)
#'
#' @return cell type marker list contains 5 genes version and 15 genes version.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindAllMarkers
#'
#' @export
#'
#'
diffmarker = function(train,
                      label_train,
                      only.pos = F){

  reference = CreateSeuratObject(train)
  reference = NormalizeData(reference)
  reference@active.ident = as.factor(label_train)
  markers = FindAllMarkers(reference,only.pos = only.pos)
  celltype = unique(markers$cluster)
  marker_5 = lapply(celltype,function(x){
    markers[markers[,6] == x,7][1:5]
  })
  marker_15 = lapply(celltype,function(x){
    markers[markers[,6] == x,7][1:15]
  })
  names(marker_5) = celltype
  names(marker_15) = celltype
  result = list()
  result$marker_5 = marker_5
  result$marker_15 = marker_15

  return(result)

}


#' marker list to matrix
#'
#' @param markers A list contains the marker gene of each cell type, the list name is cell type.
#'
#' @return A 0-1 matrix that row is marker gene, col is cell type.
#'
#' @export
#'
#'
markerconvert = function(markers){

  genes = sort(unique(unlist(markers)))
  marker_matrix = as.data.frame(matrix(0,length(genes),length(markers)))
  rownames(marker_matrix) = genes
  colnames(marker_matrix) = names(markers)
  for(i in 1:length(markers)){
    marker_matrix[markers[[i]],i] = 1
  }

  return(marker_matrix)

}

#' marker list to digitalcellsorter matrix
#'
#' @param markers A list contains the marker gene of each cell type, the list name is cell type.
#'
#' @return A 0-1 matrix that row is marker gene, col is cell type.
#'
#' @export
#'
#'
digitalcellsorter_matrix = function(markers){

  marker_matrix = markerconvert(markers)
  d_matrix = rbind(Marker=colnames(marker_matrix),marker_matrix,stringsAsFactors=F)

  return(d_matrix)

}


