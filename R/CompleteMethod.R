#'------------------------
#'---Supervised methods---
#'------------------------
#'
#' ACTINN method
#' The tensorflow1 is used in ACTINN, please set the right python version in R through reticulate::use_python().
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param python_link The absolute path of python that tensorflow 1 is installed.
#'
#' @return The predict cell type.
#'
#' @importFrom rhdf5 h5createFile
#'
#'
#' @export
#'
#'
actinn = function(train,
                  test,
                  label_train,
                  python_link = NULL){

  data_dir = system.file('extdata/ACTINN',package='RoboticCellTypeIdentification')

  train_csv = file.path(data_dir,'train.csv')
  test_csv = file.path(data_dir,'test.csv')
  label_file = file.path(data_dir,'label.csv')
  write.table(train,train_csv,quote=F,sep=',')
  write.table(test,test_csv,quote=F,sep=',')
  write.table(data.frame(colnames(train),unname(label_train)),label_file,col.names = F,quote=F,row.names=F,sep='\t')

  train_h5 = file.path(data_dir,'train.h5')
  test_h5 = file.path(data_dir,'test.h5')
  readdata.py = file.path(data_dir,'actinn_format.py')
  actinn.py = file.path(data_dir,'actinn_predict.py')
  predictfile = file.path(data_dir,'predict.txt')

  if(is.null(python_link)){
    system(paste0('python ',readdata.py,' -i ',train_csv,' -o ',train_h5,' -f csv'))
    system(paste0('python ',readdata.py,' -i ',test_csv,' -o ',test_h5,' -f csv'))
    start_time = Sys.time()
    system(paste0('python ',actinn.py,' -trs ',train_h5,' -trl ',label_file,' -ts ',test_h5,' -o ',predictfile))
    end_time = Sys.time()
    print('Running time')
    print(as.numeric(difftime(end_time,start_time,units = 'secs')))
  }else{

    system(paste0(python_link,' ',readdata.py,' -i ',train_csv,' -o ',train_h5,' -f csv'))
    system(paste0(python_link,' ',readdata.py,' -i ',test_csv,' -o ',test_h5,' -f csv'))
    start_time = Sys.time()
    system(paste0(python_link,' ',actinn.py,' -trs ',train_h5,' -trl ',label_file,' -ts ',test_h5,' -o ',predictfile))
    end_time = Sys.time()
    print('Running time')
  }

  predict_label = read.table(predictfile,head=T,sep='\t',stringsAsFactors = F,check.names = F)[,2]

  if(file.exists(train_csv)){file.remove(train_csv)}
  if(file.exists(test_csv)){file.remove(test_csv)}
  if(file.exists(train_h5)){file.remove(train_h5)}
  if(file.exists(test_h5)){file.remove(test_h5)}
  if(file.exists(label_file)){file.remove(label_file)}
  if(file.exists(predictfile)){file.remove(predictfile)}

  return(predict_label)

}

#' CaSTLe method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param threshold threshold for 'Unknown' cell type.
#'
#' @return The predict cell type.
#'
#' @importFrom igraph compare
#' @importFrom xgboost xgboost
#'
#' @export
#'
#'
castle = function(train,
                  test,
                  label_train,
                  threshold = 0.7){

  library(igraph)

  train = t(train)
  test = t(test)

  label_train = as.factor(label_train)

  BREAKS = c(-1, 0, 1, 6, Inf)
  nFeatures = 100

  train_nfeature = colSums(train>0)
  test_nfeature = colSums(test>0)

  counts = rbind(train, test)

  is_train = c(rep(TRUE,nrow(train)), rep(FALSE,nrow(test)))

  top_mean = colnames(counts)[order(apply(counts, 2, mean), decreasing = T)]

  celltype_number = length(levels(label_train))

  targetClassification = as.data.frame(matrix(rep(0,celltype_number*sum(!is_train)), nrow=celltype_number), row.names = levels(label_train))

  for (cell_type in levels(label_train)) {

    cell_type_binary = as.factor(ifelse(label_train == cell_type, cell_type, paste0("NOT",cell_type)))

    top_mutual = names(sort(apply(train,2,function(x) { compare(cut(x,breaks=BREAKS),cell_type_binary,method = "nmi") }), decreasing = T))

    selected_feature = union(head(top_mean, nFeatures) , head(top_mutual, nFeatures) )

    tmp = cor(counts[,selected_feature], method = "pearson")
    tmp[!lower.tri(tmp)] = 0

    selected_feature = selected_feature[apply(tmp,2,function(x) any(x < 0.9))]

    dsBins = apply(counts[, selected_feature], 2, cut, breaks= BREAKS)

    nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })

    ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))

    ds0 = ds0[,-1]

    cat(paste0("<h2>Classifier for ",cell_type,"</h2>"))

    inTypeSource = label_train == cell_type

    xg = xgboost(data=ds0[is_train,] ,
                 label=inTypeSource,
                 objective="binary:logistic",
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)

    inTypeProb = predict(xg, ds0[!is_train, ])

    targetClassification[cell_type,] = inTypeProb

  }

  predict_label = rownames(targetClassification)[apply(targetClassification,2,which.max)]
  if(threshold > 0){
    predict_label[which(apply(targetClassification,2,max) < threshold)] = "Unknown"
  }

  return(predict_label)

}

#' CEHTAH method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param threshold threshold for 'Unknown' cell type.
#'
#' @return The predict cell type.
#'
#' @importFrom CHETAH CHETAHclassifier
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
#'
#'
chetah = function(train,
                  test,
                  label_train){

  sce = SingleCellExperiment(assays = list(counts = train),colData = data.frame(celltypes = label_train))
  sce_test = SingleCellExperiment(assays = list(counts = test))

  sce_test = CHETAHclassifier(input = sce_test, ref_cells = sce)

  predict_label = unname(sce_test$celltype_CHETAH)

  return(predict_label)

}

#' garnett method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param marker_list A list contains marker gene of cell type.
#' @param cds Input CDS object.
#' @param marker_file A character path to the marker file to define cell types.See details and documentation for \code{\link{Parser}} by running \code{?Parser}for more information.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs. For example, for humans use org.Hs.eg.db. See available packages at \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}. If your organism does not have an AnnotationDb-class database available, you can specify "none", however then Garnett will not check/convert gene IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if db = "none".
#' @param marker_file_gene_id_type The type of gene ID used in the marker file. Should be one of the values in \code{columns(db)}. Default is "SYMBOL". Ignored if db = "none".
#' @param min_observations An integer. The minimum number of representative cells per cell type required to include the cell type in the predictive model. Default is 8.
#' @param max_training_samples An integer. The maximum number of representative cells per cell type to be included in the model training. Decreasing this number increases speed, but may hurt performance of the model. Default is 500.
#' @param num_unknown An integer. The number of unknown type cells to use as an outgroup during classification. Default is 500.
#' @param propogate_markers Logical. Should markers from child nodes of a cell type be used in finding representatives of the parent type? Should generally be \code{TRUE}.
#' @param cores An integer. The number of cores to use for computation.
#' @param lambdas \code{NULL} or a numeric vector. Allows the user to pass their own lambda values to \code{\link[glmnet]{cv.glmnet}}. If \code{NULL}, preset lambda values are used.
#' @param classifier_gene_id_type The type of gene ID that will be used in the classifier. If possible for your organism, this should be "ENSEMBL", which is the default. Ignored if db = "none".
#' @param return_initial_assign Logical indicating whether an initial assignment data frame for the root level should be returned instead of a classifier. This can be useful while choosing/debugging markers. Please note that this means that a classifier will not be built, so you will not be able to move on to the next steps of the workflow until you rerun the functionwith \code{return_initial_assign = FALSE}. Default is \code{FALSE}.
#' @param classifier Trained garnett_classifier - output from \code{\link{train_cell_classifier}}.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if db = "none".
#' @param rank_prob_ratio Numeric value greater than 1. This is the minimum odds ratio between the probability of the most likely cell type to the second most likely cell type to allow assignment. Default is 1.5. Higher values are more conservative.
#' @param cluster_extend Logical. When \code{TRUE}, the classifier provides a secondary cluster-extended classification, which assigns type for the entire cluster based on the assignments of the cluster members. If the pData table of the input CDS has a column called "garnett_cluster", this will be used for cluster-extended assignments. Otherwise, assignments are calculated using Louvain community detection in PCA space. This assignment is returned as a column in the output CDS pData table. For large datasets, if the "garnett_cluster" column is not provided and \code{cluster_extend = TRUE}, the function can be significantly slower the first time it is run. See details for more information.
#' @param verbose Logical. Should progress messages be printed.
#' @param cluster_extend_max_frac_unknown Numeric between 0 and 1. The maximum fraction of a cluster allowed to be classified as 'Unknown' and still extend classifications to the cluster. Only used when \code{cluster_extend = TRUE}. Default is 0.95. See details.
#' @param cluster_extend_max_frac_incorrect Numeric between 0 and 1. The maximum fraction of classified cells in a cluster allowed to be incorrectly classified (i.e. assigned to a non-dominant type) and still extend classifications to the cluster. Fraction does not include 'Unknown' cells. Only used when \code{cluster_extend = TRUE}. Default is 0.1. See details.
#' @param return_type_levels Logical. When \code{TRUE}, the function additionally appends assignments from each hierarchical level in the classifier as columns in the pData table labeled \code{cell_type_li}, where "i" indicates the corresponding level index
#'
#'
#' @return The predict cell type.
#'
#' @importFrom monocle newCellDataSet
#' @importFrom garnett train_cell_classifier classify_cells
#' @importFrom BiocGenerics estimateSizeFactors
#' @importFrom Biobase pData
#' @import org.Mm.eg.db
#'
#'
#' @export
#'
#'
garnett = function(train,
                   test,
                   label_train,
                   marker_list,
                   db = org.Mm.eg.db,
                   cds_gene_id_type = "SYMBOL",
                   marker_file_gene_id_type = "SYMBOL",
                   min_observations = 8,
                   max_training_samples = 500,
                   num_unknown = 500,
                   propogate_markers = T,
                   cores = 1,
                   lambdas = NULL,
                   return_initial_assign = F,
                   rank_prob_ratio = 1.5,
                   cluster_extend = FALSE,
                   verbose = FALSE,
                   cluster_extend_max_frac_unknown = 0.95,
                   cluster_extend_max_frac_incorrect = 0.1,
                   return_type_levels = FALSE){

  fdata = data.frame(gene_short_name = rownames(train),num_cells_expressed = rowSums(train>0))
  rownames(fdata) = rownames(train)
  pdata = data.frame(FACS_type = label_train)
  rownames(pdata) = colnames(train)

  pd = new("AnnotatedDataFrame", data = pdata)
  fd = new("AnnotatedDataFrame", data = fdata)
  cds_train = newCellDataSet(as(train, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)

  cds_train = estimateSizeFactors(cds_train)

  marker_dir = system.file('extdata',package='RoboticCellTypeIdentification')
  marker_file = file.path(marker_dir,'marker.txt')

  for(i in 1:length(marker_list)){
    output1 = paste0('>',names(marker_list)[i])
    output2 = paste0('expressed: ',paste(marker_list[[i]],collapse = ', '))
    write(output1,file=marker_file,append = T)
    write(output2,file=marker_file,append = T)
    write('',file=marker_file,append = T)
  }

  garnett_classifier = train_cell_classifier(cds = cds_train,
                                             marker_file = marker_file,
                                             db = db,
                                             cds_gene_id_type = cds_gene_id_type,
                                             marker_file_gene_id_type = marker_file_gene_id_type,
                                             min_observations = min_observations,
                                             max_training_samples = max_training_samples,
                                             num_unknown = num_unknown,
                                             propogate_markers = propogate_markers,
                                             cores = cores,
                                             lambdas = lambdas,
                                             return_initial_assign = return_initial_assign)

  file.remove(marker_file)

  fdata = data.frame(gene_short_name = rownames(test),num_cells_expressed = rowSums(test>0))
  rownames(fdata) = rownames(test)
  fd = new("AnnotatedDataFrame", data = fdata)
  cds_test = newCellDataSet(test,featureData = fd)
  cds_test = estimateSizeFactors(cds_test)

  cds_test = classify_cells(cds = cds_test,
                            classifier = garnett_classifier,
                            db = db,
                            cds_gene_id_type = cds_gene_id_type,
                            rank_prob_ratio = rank_prob_ratio,
                            cluster_extend = cluster_extend,
                            verbose = verbose,
                            cluster_extend_max_frac_unknown = cluster_extend_max_frac_unknown,
                            cluster_extend_max_frac_incorrect = cluster_extend_max_frac_incorrect,
                            return_type_levels = return_type_levels)

  predict_label = pData(cds_test)$cell_type

  return(predict_label)

}

#' scibet method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param norm The data need normalize or not.
#'
#' @return The predict cell type.
#'
#' @importFrom scibetR SciBet_R
#'
#' @export
#'
#'
scibet = function(train,
                  test,
                  label_train,
                  norm = F){

  label_train = as.factor(label_train)

  if(!norm){
    train = log2(train + 1)
    test = log2(test + 1)
  }

  exp = cbind(as.data.frame(t(train)),label_train)

  colnames(exp)[ncol(exp)] = 'label'

  test = t(test)

  predict_label = SciBet_R(exp,test)

  return(predict_label)

}

#' scID method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param model Classification model supported via ‘caret’ package. A list of all models can be found here.
#' @param reclassify Cell types to reclassify using a different model
#'
#' @return The predict cell type.
#'
#' @importFrom scID scid_multiclass
#'
#' @export
#'
#'
scid = function(train,
                test,
                label_train,
                logFC = 0.6,
                only_pos = FALSE,
                estimate_weights_from_target = FALSE){

  scID_output = scid_multiclass(target_gem = test,reference_gem = train,reference_clusters = label_train,logFC = logFC,only_pos = only_pos,estimate_weights_from_target = estimate_weights_from_target)

  predict_label = unname(scID_output$labels)

  return(predict_label)

}

#' scLearn method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param species The species the cells belong to. Currently species "Hs" for homo sapiens or species "Mm" for mus musculus are available. It is used to detect mitochondrial genes, so species doesn't matter if data won't be considered the percentage of mitochondrial genes.
#' @param gene_low The minimum gene number for cells. Cells with genes below this threshold are filtered(default 500).
#' @param gene_high The maximum gene number for cells. cells with genes above this threshold are filtered(default 10000).
#' @param mito_high The maximum percentage of mitochondrial genes detected(default 0.1 among all the detected genes).
#' @param umi_low The minimum number of unique molecular identifier induced (default 1500).
#' @param umi_high The maximum number of unique molecular identifier induced (default Inf).
#' @param logNormalize TRUE by default. If FALSE, the data will not be performed log normalized.
#' @param log_normalized Boolean, whether the train is log normalized (default : False).
#' @param threshold threshold to select feature, detail in M3Drop package.
#' @param sample_information_timePoint A character vector showing the time point of each sample. The column name of the vector is the sample name (default:NULL).
#' @param bootstrap_times The times for bootstrapping which should be at least larger than 1 (default:10).
#' @param cutoff The cutoff for selecting similarity threshold for each cell type (default:0.01).
#' @param dim_para The threshold to choose proper dimension for MDDM (default:0.999).
#' @param vote_rate Default is 0.6.
#' @param diff Default is 0.05.
#' @param threshold_use Default is TRUE.
#'
#' @return The predict cell type.
#'
#' @importFrom scLearn Feature_selection_M3Drop scLearn_model_learning scLearn_cell_assignment Cell_qc Cell_type_filter
#'
#' @export
#'
#'
sclearn = function(train,
                   test,
                   label_train,
                   species = 'Mm',
                   gene_low = 500,
                   gene_high = 10000,
                   mito_high = 0.1,
                   umi_low = 1500,
                   umi_high = Inf,
                   logNormalize = TRUE,
                   log_normalized = F,
                   min_cell_number = 10,
                   threshold = 0.05,
                   sample_information_timePoint = NULL,
                   bootstrap_times = 10,
                   cutoff = 0.01,
                   dim_para = 0.999,
                   vote_rate = 0.6,
                   diff = 0.05){

  train_qc = Cell_qc(train,label_train,species = species,gene_low = gene_low,gene_high = gene_high,mito_high = mito_high,umi_low = umi_low,umi_high = umi_high, logNormalize = logNormalize)

  data_type_filtered = Cell_type_filter(train_qc$expression_profile,train_qc$sample_information_cellType,min_cell_number = min_cell_number)

  selected_feature = Feature_selection_M3Drop(data_type_filtered$expression_profile,log_normalized = log_normalized,threshold = threshold)

  model = scLearn_model_learning(selected_feature,data_type_filtered$expression_profile,data_type_filtered$sample_information_cellType,sample_information_timePoint=sample_information_timePoint,bootstrap_times = bootstrap_times,cutoff = cutoff,dim_para = dim_para)

  predict_label = scLearn_cell_assignment(model,test,vote_rate = vote_rate, diff = diff)

  predict_label = predict_label[,2]

  return(predict_label)

}

#' scmapcluster method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param threshold threshold on similarity (or probability for SVM and RF)
#'
#' @return The predict cell type.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts logcounts<-
#' @importFrom scmap selectFeatures indexCluster scmapCluster
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
#'
#'
scmapcluster = function(train,
                        test,
                        label_train,
                        threshold = 0.7){

  sce = SingleCellExperiment(list(counts = train),colData = data.frame(cell_type1 = label_train))
  logcounts(sce) = log2(counts(sce) + 1)
  rowData(sce)$feature_symbol = rownames(sce)
  sce = selectFeatures(sce)

  sce_test = SingleCellExperiment(list(counts = test))
  logcounts(sce_test) = log2(counts(sce_test) + 1)
  rowData(sce_test)$feature_symbol = rownames(sce_test)

  sce = indexCluster(sce)
  scmapCluster_results = scmapCluster(projection = sce_test,index_list = list(sce@metadata$scmap_cluster_index),threshold = threshold)
  predict_label = scmapCluster_results$combined_labs

  return(predict_label)

}

#' scPred method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param model Classification model supported via ‘caret’ package. A list of all models can be found here.
#' @param reclassify Cell types to reclassify using a different model
#'
#' @return The predict cell type.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP
#' @importFrom scPred getFeatureSpace trainModel scPredict
#' @importFrom magrittr %>%
#'
#' @export
#'
#'
scpred = function(train,
                  test,
                  label_train,
                  model = 'svmRadial',
                  reclassify = NULL){

  reference = CreateSeuratObject(train)
  query = CreateSeuratObject(test)

  reference = reference %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)

  query = NormalizeData(query)

  reference$cell_type = label_train

  reference = getFeatureSpace(reference, "cell_type")

  reference = trainModel(reference,model = model,reclassify = reclassify)

  query = scPredict(query, reference)

  predict_label = unname(query$scpred_prediction)

  return(predict_label)

}

#' scVI method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param python_link The absolute path of python that tensorflow 1 is installed.
#'
#' @return The predict cell type.
#'
#' @importFrom rhdf5 h5createFile
#'
#'
#' @export
#'
#'
scvi = function(train,
                test,
                label_train,
                python_link = NULL,
                n_epochs = 20){


  scvi = reticulate::import('scvi')

  sc = reticulate::import('scanpy')

  data_dir = system.file('extdata',package='RoboticCellTypeIdentification')

  combine_csv = file.path(data_dir,'combine_data.csv')
  combine_label_csv = file.path(data_dir,'combine_label.csv')

  combine_data = t(cbind(train,test))

  combine_label = as.factor(unname(c(label_train,rep('test',ncol(test)))))
  code = levels(combine_label)
  combine_label = as.factor(as.integer(combine_label))
  names(code) = levels(combine_label)

  write.table(combine_data,combine_csv,quote=F,sep=',')
  write.table(combine_label,combine_label_csv,quote=F,sep=',')

  combine_scvi = scvi$dataset$CsvDataset(combine_csv,save_path='',sep=',',labels_file=combine_label_csv,gene_by_cell = F, new_n_genes = F)
  scanvi = scvi$models$SCANVI(combine_scvi$nb_genes, combine_scvi$n_batches, combine_scvi$n_labels)
  trainer_scanvi = scvi$inference$SemiSupervisedTrainer(scanvi, combine_scvi)

  train_indice = (0:(ncol(train)-1))
  test_indice = ncol(train):(ncol(train)+ncol(test)-1)
  trainer_scanvi$labelled_set = trainer_scanvi$create_posterior(indices=train_indice, shuffle = F)
  trainer_scanvi$labelled_set$to_monitor = c('ll','accuracy')
  trainer_scanvi$unlabelled_set = trainer_scanvi$create_posterior(indices=test_indice, shuffle = F)
  trainer_scanvi$unlabelled_set$to_monitor = c('ll','accuracy')

  trainer_scanvi$train(as.integer(n_epochs))

  result = trainer_scanvi$unlabelled_set$compute_predictions()

  predict_label = unname(code[as.numeric(result[[2]])+1])

  return(predict_label)

}

#' Seurat method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#'
#' @return The predict cell type.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures GetAssayData VariableFeatures
#'
#' @export
#'
#'
seurat = function(train,
                  test,
                  label_train,
                  k.filter = 200){

  reference = CreateSeuratObject(train)
  reference = NormalizeData(reference)
  reference = FindVariableFeatures(reference)
  reference$celltype = label_train
  query = CreateSeuratObject(test)
  query = NormalizeData(query)
  query = FindVariableFeatures(query)
  if(min(ncol(reference),ncol(test))<200){
    k.filter=50
    print(1)
  }
  if(min(ncol(reference),ncol(test))<50){
    return(NULL)
  }
  anchors = FindTransferAnchors(reference,query,k.filter = k.filter)
  predictions = TransferData(anchors,reference$celltype)
  query = AddMetaData(query, metadata = predictions)
  predict_label = unname(query$predicted.id)

  return(predict_label)

}

#' singleCellNet method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param colName_samp the name of the column that contains sample names (default : "row.names").
#' @param nTopGenes the number of classification genes per category (default : 10).
#' @param nTopGenePairs the number of top gene pairs per category (default : 25).
#' @param nRand number of random profiles generate for training (default : 70).
#' @param nTrees number of trees for random forest classifier (default : 1000).
#' @param nrand the number of random profiles generate for evaluation process (default : 50).
#'
#' @return The predict cell type.
#'
#' @importFrom singleCellNet scn_train scn_predict assign_cate
#'
#' @export
#'
#'
singlecellnet = function(train,
                         test,
                         label_train,
                         colName_samp="row.names",
                         nTopGenes = 10,
                         nTopGenePairs = 25,
                         nRand = 70,
                         nTrees = 1000,
                         nrand = 50){

  label_train = as.data.frame(label_train)

  class_info = scn_train(stTrain = label_train, expTrain = train, nTopGenes = nTopGenes, nRand = nRand, nTrees = nTrees, nTopGenePairs = nTopGenePairs, dLevel = "label_train", colName_samp = colName_samp)

  classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=test, nrand = nrand)

  predict_label = c()
  predict_label = unlist(assign_cate(classRes_val_all,predict_label))
  predict_label = predict_label[-c((length(predict_label)-49):length(predict_label))]

  return(predict_label)

}

#' SingleR method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param method Deprecated.
#' @param clusters A character vector or factor of cluster identities for each cell in \code{test}. If set, annotation is performed on the aggregated cluster profiles, otherwise it defaults to per-cell annotation.
#' @param genes,sd.thresh Arguments controlling the choice of marker genes used for annotation, see \code{\link{trainSingleR}}.
#' @param quantile,fine.tune,tune.thresh,prune Further arguments to pass to \code{\link{classifySingleR}}.
#' @param assay.type.test An integer scalar or string specifying the assay of \code{test} containing the relevant expression matrix, if \code{test} is a \linkS4class{SummarizedExperiment} object.
#' @param assay.type.ref An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix, if \code{ref} is a \linkS4class{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#'
#' @return The predict cell type.
#'
#' @importFrom SingleR SingleR
#' @importFrom Seurat CreateSeuratObject NormalizeData GetAssayData
#'
#' @export
#'
#'
singler = function(train,
                   test,
                   label_train,
                   method = 'single',
                   clusters = NULL,
                   genes = "de",
                   sd.thresh = 1,
                   quantile = 0.8,
                   fine.tune = TRUE,
                   tune.thresh = 0.05,
                   prune = TRUE,
                   assay.type.test = "logcounts",
                   assay.type.ref = "logcounts",
                   check.missing=TRUE){

  train_object = CreateSeuratObject(train)
  test_object = CreateSeuratObject(test)
  train_object = NormalizeData(train_object,verbose = F)
  test_object = NormalizeData(test_object,verbose = F)
  train = as.matrix(GetAssayData(train_object,slot = 'data'))
  test = as.matrix(GetAssayData(test_object,slot = 'data'))

  model = SingleR(test = test,
                  ref = train,
                  labels = label_train,
                  method = method,
                  clusters = clusters,
                  genes = genes,
                  sd.thresh = sd.thresh,
                  quantile = quantile,
                  fine.tune = fine.tune,
                  tune.thresh = tune.thresh,
                  prune = prune,
                  assay.type.test = assay.type.test,
                  assay.type.ref = assay.type.ref,
                  check.missing = check.missing)

  predict_label = model$labels

  return(predict_label)

}
#'------------------------
#'---Supervised methods---
#'------------------------


#'--------------------------
#'---UnSupervised methods---
#'--------------------------
#'
#' scmapcell method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param w1 a positive integer specifying the number of nearest neighbours to find.
#' @param w2 an integer specifying the number of nearest neighbours to find.
#' @param threshold the threshold which the maximum similarity between the query and a reference cell must exceed for the cell-type to be assigned.
#'
#' @return The predict cell type.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts logcounts
#' @importFrom scmap selectFeatures indexCell scmapCell scmapCell2Cluster
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#'
#'
#' @export
#'
#'
scmapcell = function(train,
                     test,
                     label_train,
                     w1 = 10,
                     w2 = 3,
                     threshold = 0.5){

  sce = SingleCellExperiment(list(counts = train),colData = data.frame(cell_type1 = label_train))
  logcounts(sce) = log2(counts(sce) + 1)
  rowData(sce)$feature_symbol = rownames(sce)
  sce = selectFeatures(sce)

  sce_test = SingleCellExperiment(list(counts = test))
  logcounts(sce_test) = log2(counts(sce_test) + 1)
  rowData(sce_test)$feature_symbol = rownames(sce_test)

  sce = indexCell(sce)
  scmapCell_results = scmapCell(sce_test,list(sce@metadata$scmap_cell_index),w = w1)
  scmapCell_clusters = scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)),w = w2,threshold = threshold)
  predict_label = scmapCell_clusters$combined_labs

  return(predict_label)

}

#' CELLBLAST method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#'
#' @return The predict cell type.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures GetAssayData VariableFeatures
#'
#' @export
#'
#'
cellblast = function(train,
                     test,
                     label_train,
                     python_link = NULL){

  # reticulate::use_condaenv('cellblast',require=T)

  if(!is.null(python_link)){
    reticulate::use_python(python_link)
  }

  cb = reticulate::import('Cell_BLAST')

  dataframe_to_cb = function(dataframe,lable=NULL){

    object = CreateSeuratObject(dataframe)
    object = NormalizeData(object)
    object = FindVariableFeatures(object,selection.method='mvp',binning.method='equal_frequency')
    if(!is.null(lable)){
      object$celltyperaw = lable
    }
    exprs = t(as.matrix(GetAssayData(object,slot='counts')))
    obs = object[[]]
    var = GetAssay(object)[[]]
    uns = VariableFeatures(object)
    data_obj = cb$data$ExprDataSet(exprs, obs, var, uns)
    return(data_obj)

  }

  data_obj_train = dataframe_to_cb(train,label_train)
  data_obj_test = dataframe_to_cb(test)

  start_time = Sys.time()
  models = list()
  for(i in 1:4){
    models[[i]] = cb$directi$fit_DIRECTi(
      data_obj_train, genes=data_obj_train$uns,
      latent_dim=as.integer(10), cat_dim=as.integer(20), random_seed=as.integer(i)
    )
  }

  blast = cb$blast$BLAST(models, data_obj_train)

  test_hits = blast$query(data_obj_test,n_neighbors=as.integer(10))$reconcile_models()$filter(by="pval", cutoff=0.05)
  predict_label = unname(unlist(test_hits$annotate("celltyperaw")))
  end_time = Sys.time()
  print(as.numeric(difftime(end_time,start_time,units = 'secs')))

  return(predict_label)

}

#' Cellfishing.ji method
#'
#' @param train A gene-cell matrix of single cell expression data.
#' @param test Same format like train.
#' @param label_train A vector of cell type label of train.
#' @param k a positive integer specifying the number of nearest neighbours to find.
#'
#' @return The predict cell type.
#'
#' @importFrom JuliaCall julia_setup
#'
#' @export
#'
#'
cellfishing.ji = function(train,
                          test,
                          label_train,
                          k = 10){

  julia = julia_setup()

  julia$library("CellFishing")
  julia$library("TableReader")

  julia$assign("train",train)
  julia$assign("test",test)
  julia$assign("featurenames_train",rownames(train))
  julia$assign("cellnames",colnames(train))
  julia$assign("featurenames_test",rownames(test))
  julia$assign("k",as.integer(k))

  julia$eval("@elapsed features = CellFishing.selectfeatures(train,featurenames_train)")
  julia$eval("@elapsed database = CellFishing.CellIndex(train,features,metadata=cellnames)")
  julia$eval("@elapsed neighbors = CellFishing.findneighbors(k, test, featurenames_test, database)")

  julia_result = julia$eval("neighbors.indexes")
  colnames(julia_result) = colnames(test)
  julia_result = apply(julia_result,2,function(x){label_train[x]})
  predict_label = apply(julia_result,2,function(x){
    re = unique(x)
    if(length(re) == 1){
      return(re)
    }else{
      return("unassigned")
    }
  })
  predict_label = unname(predict_label)
  return(predict_label)

}
#'--------------------------
#'---UnSupervised methods---
#'--------------------------


#'----------------------------
#'---SemiSupervised methods---
#'----------------------------
#'
#' DigitalCellSorter method
#'
#' @param test Same format like train.
#' @param python_link The absolute path of python that DigitalCellSorter is installed, default is null.
#' @param python_package_link If DigitalCellSorter is installed and can not load, you can provide the absolute path of DigitalCellSorter package, like '/usr/local/lib/python3.8/site-packages/', default is null.
#' @param marker_list cell type marker list contains marker genes.
#'
#' @return The predict cell type.
#'
#' @importFrom xlsx write.xlsx
#'
#' @export
#'
#'
digitalcellsorter = function(test,
                             python_link = NULL,
                             python_package_link = NULL,
                             marker_list = NULL){

  if(!is.null(python_link)){
    reticulate::use_python(python_link)
  }

  if(!is.null(python_package_link)){
    dcs = reticulate::import_from_path('DigitalCellSorter',path=python_package_link)
  }else{
    dcs = reticulate::import('DigitalCellSorter')
  }

  dcs_object = dcs$DigitalCellSorter()
  test = as.data.frame(test)
  dcs_object$prepare(test)
  dcs_object$geneNamesType = 'hugo'
  dcs_object$doQualityControl = F
  dcs_object$process()
  dcs_object$minimumNumberOfMarkersPerCelltype = 1

  if(!is.null(marker_list)){
    d_matrix = digitalcellsorter_matrix(marker_list)
    data_dir = system.file('extdata',package='RoboticCellTypeIdentification')
    marker_matrix_xlsx = file.path(data_dir,'marker_matrix.xlsx')
    write.xlsx(d_matrix,marker_matrix_xlsx)
    dcs_object$geneListFileName = marker_matrix_xlsx
  }

  dcs_object$annotate()
  predict_label = dcs_object$getCells()
  predict_label = unname(gsub(' #\\d+','',predict_label))

  return(predict_label)

}

#' SCINA method
#'
#' @param test Same format like train.
#' @param marker_list A list contains marker gene of cell type.
#' @param max_iter An integer > 0. Default is 100. Max iterations allowed for the EM algorithm.
#' @param convergence_n An integer > 0. Default is 10. Stop the SCINA algorithm if during the last n rounds of iterations, cell type assignment keeps steady above the convergence_rate.
#' @param convergence_rate A float between 0 and 1. Default is 0.99. Percentage of cells for which the type assignment remains stable for the last n rounds.
#' @param sensitivity_cutoff A float between 0 and 1. Default is 1. The cutoff to remove signatures whose cells types are deemed as non-existent at all in the data by the SCINA algorithm.
#' @param rm_overlap A binary value, default 1 (TRUE), denotes that shared symbols between signature lists will be removed. If 0 (FALSE) then allows different cell types to share the same identifiers.
#' @param allow_unknown A binary value, default 1 (TRUE). If 0 (FALSE) then no cell will be assigned to the 'unknown' category.
#' @param log_file A string names the record of the running status of the SCINA algorithem, default 'SCINA.log'.
#'
#' @return The predict cell type.
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom SCINA SCINA
#'
#' @export
#'
#'
scina = function(test,
                 marker_list,
                 max_iter = 100,
                 convergence_n = 10,
                 convergence_rate = 0.999,
                 sensitivity_cutoff = 0.9,
                 rm_overlap=FALSE,
                 allow_unknown=TRUE,
                 log_file='SCINA.log'){

  test = log(test+1)
  test[] = normalize.quantiles(test)

  marker_list = sapply(marker_list,function(x) unique(x[(!is.na(x)) & (x %in% rownames(test))]),simplify=F)

  result = SCINA(test,marker_list,max_iter = max_iter,convergence_n = convergence_n,convergence_rate = convergence_rate, sensitivity_cutoff = sensitivity_cutoff, rm_overlap=rm_overlap, allow_unknown=allow_unknown, log_file=log_file)

  predict_label = result$cell_labels

  return(predict_label)

}




# #' cellassign method
# #'
# #' @param test Same format like train.
# #' @param marker_list A list contains marker gene of cell type.
# #'
# #' @return The predict cell type.
# #'
# #' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors
# #' @importFrom cellassign cellassign marker_list_to_mat
# #' @importFrom scran computeSumFactors
# #'
# #' @export
# #'
# #'
# cellassign_r = function(test,
#                         marker_list,
#                         python_link = NULL){
#
#   sce_test = SingleCellExperiment(list(counts = test))
#   s = computeSumFactors(sce_test)
#
#   marker_mat = marker_list_to_mat(marker_list)
#
#   if(!is.null(python_link)){
#     reticulate::use_python(python_link)
#   }
#
#   intersect_marker = intersect(rownames(sce_test),rownames(marker_mat))
#   marker_mat = marker_mat[intersect_marker,]
#
#   fit = cellassign(exprs_obj = sce_test[intersect_marker,],
#                    marker_gene_info = marker_mat,
#                    s = sizeFactors(s),
#                    learning_rate = 1e-2,
#                    shrinkage = TRUE,
#                    verbose = FALSE)
#
#   predict_label = fit$cell_type
#
#   return(predict_label)
#
# }


#' SCSA method
#'
#' @param test Same format like train.
#' @param marker_list A list contains marker gene of cell type.
#' @param marker_file The absolute path of marker file (now just support marker info from Seurat FindAllMarkers method).
#'
#' @return The predict cell type.
#'
#'
#' @export
#'
#'
scsa = function(test = NULL,
                marker_list = NULL,
                marker_file = NULL){

  if(!is.null(marker_file)){

    data_dir = system.file('extdata/SCSA',package='RoboticCellTypeIdentification')
    scsa.py = file.path(data_dir,'SCSA.py')
    whole.db = file.path(data_dir,'whole.db')
    scsa_result = file.path(data_dir,'scsa_result')
    scsa_go_result = file.path(data_dir,'scsa_result.go')

    system(paste0('python3 ',scsa.py,' -d ',whole.db,' -i ',marker_file,' -s seurat -E -f1.5 -p 0.01 -o ',scsa_result,' -m txt'))
    predict_label = read.table(scsa_result,head=T,sep='\t',stringsAsFactors = F)

    if(file.exists(scsa_result)){file.remove(scsa_result)}
    if(file.exists(scsa_go_result)){file.remove(scsa_go_result)}

    return(predict_label)

  }else{
    if(!is.null(test)){

      object = CreateSeuratObject(test)
      object = NormalizeData(object,verbose = F)
      object = FindVariableFeatures(object,verbose = F)
      object = ScaleData(object,verbose = F)
      object = RunPCA(object,verbose = F)
      object = FindNeighbors(object,verbose = F)
      object = FindClusters(object,verbose = F)
      marker_table = FindAllMarkers(object,min.pct = 0.25,features=VariableFeatures(object))

      data_dir = system.file('extdata/SCSA',package='RoboticCellTypeIdentification')
      output = file.path(data_dir,'marker_table.csv')
      write.csv(marker_table,output,quote=F,row.names = F)

      data_dir = system.file('extdata/SCSA',package='RoboticCellTypeIdentification')
      scsa.py = file.path(data_dir,'SCSA.py')
      whole.db = file.path(data_dir,'whole.db')
      scsa_result = file.path(data_dir,'scsa_result')
      scsa_go_result = file.path(data_dir,'scsa_result.go')

      if(is.null(marker_list)){
        system(paste0('python3 ',scsa.py,' -d ',whole.db,' -i ',output,' -s seurat -E -f1.5 -p 0.01 -o ',scsa_result,' -m txt'))
      }else{
        user.table = file.path(data_dir,'user.table')
        genelist = c()
        for(k in 1:length(marker_list)){
          genelist = rbind(genelist,data.frame(names(marker_list)[k],marker_list[[k]]))
        }
        write.table(genelist,user.table,sep='\t',row.names = F,col.names = F,quote=F)
        system(paste0('python3 ',scsa.py,' -d ',whole.db,' -i ',output,' -s seurat -E -f1.5 -p 0.01 -o ',scsa_result,' -m txt -N -M ',user.table))
      }

      scsa_predict = read.table(scsa_result,head=T,sep='\t',stringsAsFactors = F)

      file.remove(output)
      file.remove(scsa_result)
      file.remove(scsa_go_result)
      file.remove(user.table)

      clusters = unique(scsa_predict[,3])
      scsa_predict_max = scsa_predict[sapply(clusters,function(x){which(scsa_predict[,3]==x)[1]}),]
      celltypeclusters = scsa_predict_max[,1]
      names(celltypeclusters) = clusters

      idents = object@active.ident
      if(length(setdiff(levels(idents),clusters)) > 0){
        unassign = rep('unassign',length(setdiff(levels(idents),clusters)))
        names(unassign) = setdiff(levels(idents),clusters)
        celltypeclusters = append(celltypeclusters,unassign)
        celltypeclusters = celltypeclusters[order(as.numeric(names(celltypeclusters)))]
      }

      predict_label = Idents(object)
      levels(predict_label) = celltypeclusters
      predict_label = unname(as.character(predict_label))
    }
    return(predict_label)
  }

}

#' Evaluate the prediction performance, source from PMID: 31500660.
#'
#' @param ture_label A list contains the ture label of cell type.
#' @param pred_label A list contains the predict label of cell type.
#'
#' @return A list contaions the table info of ture/predict label, the F1 score of each cell type, the median F1 score predict label, the mean F1 score of
#' predict label, the accuracy rate of predict label,the percentage of 'Unknown' cell type, the table info of ture label.
#'
#' @importFrom raster rowSums colSums
#' @importFrom stats median
#'
#' @export
#'
Evaluate = function(ture_label,
                    pred_label){

  ture_label = as.character(unlist(ture_label))
  pred_label = as.character(unlist(pred_label))
  unique_true = unlist(unique(ture_label))
  unique_pred = unlist(unique(pred_label))
  unique_all = unique(c(unique_true,unique_pred))
  conf = table(ture_label,pred_label)
  pop_size = rowSums(conf)
  pred_label = gsub('Node..','Node',pred_label)
  conf_F1 = table(ture_label,pred_label,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown','rejected','other','unassign'))
  F1 = vector()
  sum_acc = 0

  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec = conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec = conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] = (2*prec*rec) / (prec + rec)
      }
      sum_acc = sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }

  pop_size = pop_size[pop_size > 0]
  names(F1) = names(pop_size)
  med_F1 = median(F1)
  mean_F1 = mean(F1)
  total = length(pred_label)
  num_unlab = sum(pred_label == 'unassigned') + sum(pred_label == 'Unassigned') + sum(pred_label == 'rand') + sum(pred_label == 'Unknown') + sum(pred_label == 'unknown') + sum(pred_label == 'Node') + sum(pred_label == 'ambiguous') + sum(pred_label == 'rejected') + sum(pred_label == 'other') + sum(pred_label == 'unassign')
  per_unlab = num_unlab / total
  acc = sum_acc/sum(conf_F1)

  result = list(Conf = conf, Med_F1 = med_F1,Mean_F1 = mean_F1, F1 = F1, Acc = acc, PercUnl = per_unlab, PopSize = pop_size)

  return(result)
}





#' Self projection of Supervised methods and Unsupervised methods
#'
#' @param dataset A gene-cell matrix to do 5 fold cross validation.
#' @param label_train A vector of cell type label of train.
#' @param method Choosing the Supervised/Unsupervised methods to assign cell type.
#' @param n Split vector to n-fold (deufault : 5).
#' @param ... other paramter used in method/classifer, the detial could see the corresponding method.
#'
#' @return The prediction of each methods.
#'
#' @export
#'
#' @importFrom caret createFolds
#'
selfprojection = function(dataset,
                          label_dataset,
                          method,
                          n = 5,
                          ...){

  foldlist = createFolds(label_dataset,n)

  evaluation_result = list()

  for(i in 1:length(foldlist)){

    train = dataset[,unlist(foldlist[-i],use.names = F)]
    test = dataset[,foldlist[[i]]]

    label_train = label_dataset[unlist(foldlist[-i],use.names = F)]
    label_test = label_dataset[foldlist[[i]]]

    method = tolower(method)
    if(method %in% c('cellfishing.ji','scmapcell','cellblast')){
      predict_label = unsupervised(train,test,label_train,method,...)
    }else{
      predict_label = supervised(train,test,label_train,method,...)
    }

    evaluation_result[[i]] = Evaluate(label_test,predict_label)

  }

  return(evaluation_result)

}







