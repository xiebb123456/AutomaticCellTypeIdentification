import scipy.stats
from scipy.spatial import distance_matrix
import math, copy
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import cluster, mixture

STR_UNASSIGNED = 'unassigned'
MIN_OFREQ = 1e-20

def get_occurrence_freq(dfX, cell_type_X, cell_type_lst_to_test = []):
    
    genes = dfX.columns.values
    N_cells = dfX.shape[0]
    b_csel = np.full( N_cells, False )
    
    if len(cell_type_lst_to_test) == 0:
        # b = ~cell_type_X.isnull()
        b = ~pd.isnull(cell_type_X)
        cell_type_lst_to_test = list(set(cell_type_X[b]))
    
    for c in cell_type_lst_to_test:
        b = cell_type_X == c
        b_csel = b_csel | b
        
    X = dfX.loc[b_csel,:]
    cell_type = np.array(cell_type_X)[b_csel]
    
    df = pd.DataFrame(index = cell_type_lst_to_test, columns = genes)
    bx = []
    for c in cell_type_lst_to_test:
        b = cell_type == c
        if np.sum(b) == 0:
            bx.append(False)
        else:
            bx.append(True)
        # X_sel = (X.loc[b,:] > 0)
        ofreq = (X.loc[b,:] > 0).sum(axis = 0).astype(np.float32)/np.sum(b)
        df.loc[c,:] = ofreq

    df_ofreq = df.loc[bx,:].copy(deep=True)
    
    return df_ofreq


def get_occurrence_freq_per_sample(dfX, cell_type_X, sample_name, cell_type_lst_to_test = []):
    
    if len(cell_type_lst_to_test) == 0:
        # b = ~cell_type_X.isnull()
        b = ~pd.isnull(cell_type_X)
        cell_type_lst_to_test = list(set(cell_type_X[b]))
    
    sample_lst = list(set(sample_name))
    sample_name_ary = np.array(sample_name)
    cnt = 0
    for sample in sample_lst:
        b = sample_name_ary == sample
        X = dfX.loc[b,:]
        cell_type = cell_type_X[b]
        df_tmp = get_occurrence_freq(X, cell_type, cell_type_lst_to_test)
        if cnt == 0:
            df_ofreq = df_tmp
            cell_type_out = list(df_tmp.index.values)
        else:
            df_ofreq = pd.concat( [df_ofreq, df_tmp], axis = 0)
            cell_type_out = cell_type_out + list(df_tmp.index.values)
        cnt += 1
        
    print(df_ofreq.shape, len(cell_type_out))
    df_ofreq = get_occurrence_freq(df_ofreq, df_ofreq.index.values, cell_type_lst_to_test)
    
    return df_ofreq 


def get_marker_profile( df_ofreq, of_th = 0.9, N_mkr = 0, stype = 'max'):
    
    genes = df_ofreq.columns.values
    ct = df_ofreq.index.values
    
    df_profile_b = df_ofreq.copy(deep=True)

    df_stat = pd.DataFrame(index = np.arange(N_mkr))
    for c in list(ct):
        b = ct == c
        of1 = df_ofreq.loc[c,:].astype(np.float32) + MIN_OFREQ
        if stype == 'max':
            of2 = df_ofreq.loc[~b,:].max(axis=0).astype(np.float32) + MIN_OFREQ
        else:
            of2 = df_ofreq.loc[~b,:].mean(axis=0) + MIN_OFREQ
        
        # score = np.log2(of1/of2)
        b1 = (of1 >= of_th).astype(int)
        # score = score*of1*b1
        score = (of1)*((1-of2)**2)
        
        if N_mkr == 0:
            # bs1 = ((score > np.log2(fold_th) & (score > 0))).astype(int)
            bf1 = (of1 >= of_th) & (of2 < (1-of_th))
        else:
            # bf1 = score > np.log2(fold_th)
            bf1 = (of1 >= of_th) & (of2 < (1-of_th)) # score > np.log2(fold_th)
            s = np.array(score)
            b = np.full(len(s), False)
            N_sel = min(N_mkr, np.sum(score > 0))
            odr = (-s).argsort()
            b[odr[:N_sel]] = True
            bs1 = b | bf1
                        
        wgt = bs1.astype(np.float32) # (of1-of2)*bs1.astype(np.float32) # 
        df_profile_b.loc[c,:] = wgt # bs1.astype(np.float32)
        
        of1_sel = of1[b]
        of2_sel = of2[b]
        score_sel = score[b]
        gene_sel = df_profile_b.columns.values[b]
        odr = (-np.array(score_sel)).argsort()
        
        df_stat[c+'_m'] = gene_sel[odr]
        df_stat[c+'_f1'] = list(of1_sel[odr]) # list(score_sel[odr])
        df_stat[c+'_f2'] = list(of2_sel[odr])
        
    return df_profile_b


def get_score( X_in, df_mkr ):
    
    cell_types = df_mkr.index.values
    X = (X_in > 0)
    genes = X.columns.values
    
    df_score = pd.DataFrame(index = X.index.values, columns = cell_types)
    for ct in cell_types:
        b = df_mkr.loc[ct,:] > 0
        gsel = df_mkr.columns.values[b]
        gsel = set(gsel).intersection(genes)
        wgt = df_mkr.loc[ct,gsel]
        # df_score.loc[:,ct] = (X.loc[:,gsel]).sum(axis=1)/np.sum(b)
        df_score.loc[:,ct] = (X.loc[:,gsel]*wgt).sum(axis=1)/np.sum(wgt)
    
    n_mkrs = df_mkr.sum(axis = 1).T
    wgt_mkrs = 1/(1+3*np.exp(-(n_mkrs-1)*0.5))    
    df_score = df_score*wgt_mkrs
            
    idx = np.array(df_score).argmax(axis=1)
    best_score = df_score.max(axis=1)
    
    cell_type_pred = cell_types[idx]
    
    return df_score, cell_type_pred, best_score


def gen_marker_profile(X, cell_type, N_mkrs, cells_to_excl = [], \
                       of_th = 0.5, ttype = 's', stype = 'mean'):
    
    df_ofreq = get_occurrence_freq(X, cell_type) 
    # df_ofreq = get_occurrence_freq_per_sample(X, cell_type, sample_name, cells_to_test)
    df_pb = get_marker_profile(df_ofreq, of_th = of_th, N_mkr = N_mkrs, stype = stype ) 

    cells_to_test = list(set(cell_type))    
    if len(cells_to_excl) > 0:
        for ct in cells_to_excl:
            if ct in cells_to_test: 
                cells_to_test.remove(ct)
        
    cells_cur = df_pb.index.values
    cells_to_test = set(cells_to_test).intersection(cells_cur) 
    b = df_pb.sum(axis = 0) > 0
    df_p = df_pb.loc[cells_to_test,b]
        
    df_s, cell_type_pred, score = get_score( X, df_p )
    
    return df_p, df_s, cell_type_pred, df_ofreq


def get_gmm_param( df_score, cell_type_true = None ):
    
    target_cell_types = list(set(df_score.columns.values))
    cells = df_score.index.values
    
    df_param = pd.DataFrame( index = target_cell_types, \
                          columns = ['w0', 'm0', 'v0', 'w1', 'm1', 'v1'] )
    b_ts = np.full(df_score.shape[0], False)
    for tct in target_cell_types:
        
        x = df_score[tct]
        
        if cell_type_true is None:
            gmm = mixture.GaussianMixture(n_components = 2, random_state = 0)
            y = gmm.fit_predict(np.array(x).reshape(-1, 1))

            mns = [m[0] for m in gmm.means_]
            cvs = [cv[0,0] for cv in gmm.covariances_]

            wgs = gmm.weights_           
            if mns[0] < mns[1]:
                w0, w1 = wgs[0], wgs[1]
                m0, m1 = mns[0], mns[1]
                v0, v1 = cvs[0], cvs[1]
            else:
                w0, w1 = wgs[1], wgs[0]
                m0, m1 = mns[1], mns[0]
                v0, v1 = cvs[1], cvs[0]
            
        else:
            b = cell_type_true == tct
            if np.sum(b) > 0:
                m1 = x[b].mean()
                v1 = x[b].var()
                w1 = np.sum(b)/len(b)
            else:
                m1 = 1
                v1 = 1e-6
                w1 = 0
                
            m0 = x[~b].mean()
            v0 = x[~b].var()
            w0 = np.sum(~b)/len(b)
            
        param = [w0, m0, v0, w1, m1, v1]
        df_param.loc[tct, :] = param
            
    return df_param
        
        
def get_normal_pdf( x, mu, var, nbins):
    
    y = np.array(x)
    mn_x = y.min()
    mx_x = y.max()
    L = 100
    # dx = len(y)*(mx_x-mn_x)/L
    dx = (mx_x-mn_x)/nbins
    xs = np.arange(mn_x,mx_x, dx )
    pdf = (dx*len(y))*np.exp(-((xs-mu)**2)/(2*var+1e-10))/(np.sqrt(2*math.pi*var)+1e-10) + 1e-10
    return pdf, xs


def get_threshold_using_param( df_param, Target_FPR = 0.1 ):
    
    target_cell_types = list(set(df_param.index.values))
    
    df_th = pd.DataFrame( index = target_cell_types, columns = ['threshold'])
    for tct in target_cell_types:
        
        w0, m0, v0, w1, m1, v1 = df_param.loc[tct,:] 
        mns = [m0, m1]
        cvs = [v0, v1]
        wgs = [w0, w1]

        pdfs = []
        for k in range(len(mns)):
            x = np.arange(1000)/1000
            pdf, xs = get_normal_pdf( x, mns[k], cvs[k], 1000 )
            pdfs.append(pdf)

        pdf0 = pdfs[0]*w0
        pdf1 = pdfs[1]*(w1 + 1e-10)

        thresh = -1
        for k, a in enumerate(xs):
            if (a > m0): # & (a < m1):
                p0 = np.sum(pdf0[k:])
                p1 = np.sum(pdf1[k:])
                if (p0/(p0+p1)) < Target_FPR:
                    thresh = a
                    break
        if thresh < 0:
            thresh = (mns[0] + mns[1])*0.5
        df_th.loc[tct,'threshold'] = thresh
            
    return df_th
                

def get_threshold( df_score, cell_type_pred, cell_type_true, \
                   Target_FPR = 0.1, verbose = False ):
    
    target_cell_types = list(set(df_score.columns.values))
    cells = df_score.index.values
    
    if verbose: 
        print('Getting threshold .. ')
    df_th = pd.DataFrame( index = target_cell_types, columns = ['threshold'])
    b_ts = np.full(len(cell_type_true), False)
    for tct in target_cell_types:
        
        score = df_score[tct]
        bt = cell_type_true == tct
        cells_t = cells[bt]
        score_t = score[bt]
        nt = np.sum(bt)
        
        b_ts = b_ts | bt
        
        bp = cell_type_pred == tct
        cells_p = cells[bp]
        score_p = score[bp]
        sp = np.mean(score_p)
        npp = np.sum(bp)
        
        # b_ts = b_ts | bp
        
        bfp = (~bt) & bp
        score_fp = score[bfp]
        sfp = np.mean(score_fp)
        nfp = np.sum(bfp)
        
        ntp = np.sum(bt & bp)
        acc = ntp/nt
        fpr = nfp/npp
        
        if nfp > (npp*Target_FPR):
            sfp_lst = list(score_fp)
            sfp_lst.sort(reverse = True)
            n = int(npp*Target_FPR)
            sth = sfp_lst[n]
            df_th.loc[tct,'threshold'] = sth
            if verbose: 
                print('%20s: ACC = %4.3f (%i/%i), FPR = %4.3f (%i/%i) -> Th = %4.3f (%4.3f, %4.3f)' \
                       % (tct, acc, ntp, nt, fpr, nfp, npp, sth, sp, sfp) )
        else:
            df_th.loc[tct,'threshold'] = score_p.min()
            if verbose: 
                print('%20s: ACC = %4.3f (%i/%i), FPR = %4.3f (%i/%i)' \
                       % (tct, acc, ntp, nt, fpr, nfp, npp) )
    
    
    cell_type_true_ts = cell_type_true[b_ts]
    cell_type_pred_ts = cell_type_pred[b_ts]
    
    if verbose:
        ## Confusion matrix & performance report
        c_mat = met.confusion_matrix(cell_type_true_ts, cell_type_pred_ts ) 
        print("Confusion Matrix: ")
        print(c_mat)
        print("Classification Report: ")
        print(met.classification_report(cell_type_true_ts, cell_type_pred_ts )) 
    
    return df_th
        
    
def correct_pred( cell_type_pred_in, best_score, df_thresh ):
    
    cell_type_pred = copy.deepcopy(cell_type_pred_in)
    cell_types = list(df_thresh.index.values)
    
    for ct in cell_types:
        b1 = cell_type_pred == ct
        b2 = best_score < df_thresh.loc[ct,'threshold']
        b = b1 & b2
        if np.sum(b) > 0:
            cell_type_pred[b] = STR_UNASSIGNED
    
    return cell_type_pred


## PCA and Clustering
def pca_and_clustering(X, N_clusters, clust_algo = 'gm', N_pca_components = 15):

    ## Get X_pca
    pca = PCA(n_components = N_pca_components, copy = True, random_state = 0)
    Xi = X #(X > 0).astype(np.float32)
    res = pca.fit(Xi)
    X_pca = Xi.dot(res.components_.T) 

    if clust_algo[:2] == 'gm':
        gmm = mixture.GaussianMixture(n_components = N_clusters, random_state = 0)
        cluster_label = gmm.fit_predict(X_pca)
    elif clust_algo[:2] == 'km':
        km = cluster.KMeans(n_clusters = N_clusters, random_state = 0)
        km.fit(X_pca)
        cluster_label = km.labels_
    else:
        cluster_label = None
        
    return cluster_label, X_pca

def print_clustering_info(cluster_label):
    
    clust_lst = list(set(cluster_label))
    n_cells = np.zeros(len(clust_lst))
    for k, c in enumerate(clust_lst):
        b = cluster_label == c
        n_cells[k] = np.sum(b)
    print('N_cells_per_cluster: max, median, min = %i, %i, %i with %i(%i) clusters(cells)' % \
          (np.max(n_cells), np.median(n_cells), np.min(n_cells), len(clust_lst), len(cluster_label)))
      

## Use GMM on the number of expressed markers and 
## use centroid based distance to assign confidence level.
def norm_mahal(row_vecs, icov): 
    row_vecs_C = row_vecs.dot(icov)
    dist_col_vec = np.sum(row_vecs_C * row_vecs, axis = 1)
    return dist_col_vec


def get_type_num( items, verbose = True ):
    
    item_lst = list(set(items))
    n_types = []
    for item in item_lst:
        b = np.array(items) == item
        n_types.append(np.sum(b))
        
    n_types = np.array(n_types)
    odr = (-n_types).argsort()
    
    if verbose:
        for o in odr:
            print('%20s: %i' % (item_lst[o], n_types[o]) )
        print('%20s: %i' % ('Median', np.median(n_types)) )
        
    item_lst = np.array(item_lst)
    return item_lst[odr], n_types[odr] 
    
def get_maj(ivec, cto, p_cells_dict, p_min = 0.1):

    items = list(set(ivec))
    if len(items) == 1:
        return cto
    
    Num = np.zeros(len(items))
    Score = np.zeros(len(items))
    for k, item in enumerate(items):
        b = ivec == item
        Num[k] = np.sum(b)
    k = np.argmax(Num)

    b = False
    if items[k] == STR_UNASSIGNED:
        odr = (-Num).argsort()        
        if len(odr) > 1:
            if Num[odr[1]] >= round(np.sum(Num)*(p_min)):
                k = odr[1]
            # elif Num[k] <= round(np.sum(Num)*(1-p_min)):
            #     return cto
       
    return  items[k]


def apply_knn(X_pca, cell_type, p_cells_dict, p_min, Nsel):

    ## Apply KNN(k-nearest neighbor) rule
    cells = X_pca.index.values
    dist_mat = distance_matrix( x = X_pca,  y = X_pca )
    for k in range( dist_mat.shape[0] ):
        dist_mat[k,k] = 1e20
     
    cell_type_new = copy.deepcopy(cell_type)
    for k in range( dist_mat.shape[0] ):
        dx = np.array(dist_mat[k,:])
        odr = dx.argsort()[:min(Nsel, len(dx))]
        ct = get_maj(cell_type[odr], cell_type[k], p_cells_dict, p_min = p_min)
        cell_type_new[k] = ct
   
    return cell_type_new
    

def revise_unassinged( cluster_label, cell_type_pred, X_pca, df_score, \
                       p_maj = 0.2, p_min = 0.2, d_red = 0.9, Nsel = 25 ):

    cluster_labels = list(set(cluster_label))    
    cell_type_pred_ary = np.array(cell_type_pred)
    cell_name_lst = list(set(cell_type_pred))
    for c in cluster_labels:
        b = cluster_label == c

        ct_sel = cell_type_pred_ary[b]
        bua = cell_type_pred_ary == STR_UNASSIGNED
        cell_name_lst, n_cells = get_type_num( ct_sel, verbose = False )
        bn = (n_cells >= 10) & (n_cells/sum(n_cells) >= 0.1 )
        cell_name_lst = cell_name_lst[bn]
        n_cells = n_cells[bn]
        p_cells_dict = dict(zip(cell_name_lst, n_cells/np.sum(n_cells)))

        if (len(n_cells) <= 1): # only one cell type
            pass
        else:
            bx = b
            cell_type_pred_ary[bx] = apply_knn(X_pca.loc[bx,:], cell_type_pred_ary[bx], p_cells_dict, p_min, Nsel)

            ## Recompute variables
            ct_sel = cell_type_pred_ary[b]               
            bua = cell_type_pred_ary == STR_UNASSIGNED
            cell_name_lst, n_cells = get_type_num( ct_sel, verbose = False )
            bn = (n_cells >= 10) & (n_cells/sum(n_cells) >= 0.1 )
            cell_name_lst = cell_name_lst[bn]
            n_cells = n_cells[bn]

            if (np.sum(b&bua) > np.sum(b)*(1-p_maj)):
                pass
            elif (np.sum(b&bua) <= np.sum(b)*(1-p_maj)): 
                
                cell_name_lst_sel = list(cell_name_lst)
                if STR_UNASSIGNED in cell_name_lst_sel:
                    cell_name_lst_sel.remove(STR_UNASSIGNED)

                if len(cell_name_lst_sel) == 1:
                    
                    m = 0
                    ct = cell_name_lst_sel[0]
                    bc = cell_type_pred_ary == ct
                    cent = X_pca.loc[b&bc,:].mean(axis = 0)
                    X_cent = X_pca.loc[b&bc,:] - cent
                    Cv = X_cent.T.dot(X_cent)/X_cent.shape[0]
                    s = np.mean(np.diagonal(Cv))
                    Cv = Cv + (s*0.1)*np.identity(Cv.shape[0])
                    icov = np.linalg.inv( Cv )
                    # sz = np.sum(b&bc)/np.sum(b)
                    cells_target = X_pca.index.values[b&bc]

                    bx = b 
                    cells_cluster = X_pca.index.values[bx]
                    X_pca_ua = X_pca.loc[bx,:]
                    cell_type_new = copy.deepcopy(cell_type_pred_ary[bx] )
                    
                    df_dist = pd.DataFrame(index = X_pca_ua.index.values, columns = cell_name_lst_sel)
                    df_dist.loc[:,:] = 1e10 
                    row_vecs = X_pca_ua - cent
                    df_dist[ct] = norm_mahal(row_vecs, icov)
                    
                    max_dist_intra = df_dist.loc[cells_target, ct].max()
                    bt = df_dist[ct] <= (max_dist_intra)*0.8
                    cell_type_new[bt] = ct
                    cell_type_pred_ary[bx] = cell_type_new                   
                    
                elif len(cell_name_lst_sel) > 1:
                    
                    centroids = []
                    iCovs = []
                    Sizes = []
                    for m, ct in enumerate(cell_name_lst_sel):
                        bc = cell_type_pred_ary == ct
                        cent = X_pca.loc[b&bc,:].mean(axis = 0)
                        X_cent = X_pca.loc[b&bc,:] - cent
                        Cv = X_cent.T.dot(X_cent)/X_cent.shape[0]
                        s = np.mean(np.diagonal(Cv))
                        Cv = Cv + (s*0.1)*np.identity(Cv.shape[0])
                        icov = np.linalg.inv( Cv )
                        sz = np.sum(b&bc)/np.sum(b)
                        centroids.append(cent)
                        iCovs.append(icov)
                        Sizes.append(sz)

                    bx = b #&bua 
                    cells = X_pca.index.values[bx]
                    X_pca_ua = X_pca.loc[bx,:]
                    df_score_ua = df_score.loc[bx,:]
                    df_dist = pd.DataFrame(index = X_pca_ua.index.values, columns = cell_name_lst_sel)
                    df_dist.loc[:,:] = 1e10 

                    for m, ct in enumerate(cell_name_lst_sel):
                        row_vecs = X_pca_ua - centroids[m]
                        df_dist[ct] = norm_mahal(row_vecs, iCovs[m]) #*(1/Sizes[m])

                    idx = np.array(df_dist).argmin(axis=1)
                    cell_type_new = df_dist.columns.values[idx]                
                    cell_type_pred_ary[bx] = cell_type_new
                
    return cell_type_pred_ary    


def print_n_mkrs( df_p, N_mkrs ):
    
    N_ms = (df_p > 0).sum(axis=1)
    s = 'N_markers: '
    for k, n in enumerate(N_ms):
        if n == N_mkrs:
            s = s + '%i(%s), ' % (n, df_p.index.values[k])
        else:
            s = s + '%i(%s), ' % (n, df_p.index.values[k])
    print(s)


def save_markers( df_mkr, file = 'sc_makers.csv', df_ofreq = None, ct_exc = [] ):
    
    cell_type_lst = list(set(df_mkr.index.values))
    
    # if 'Unknown' in cell_type_lst: cell_type_lst.remove('Unknown')
    # if 'Tumor' in cell_type_lst: cell_type_lst.remove('Tumor')
    if len(ct_exc) > 0:
        for ct in ct_exc:
            if ct in cell_type_lst: 
                cell_type_lst.remove(ct)
    
    df_mkr = df_mkr.loc[cell_type_lst,:]
    if df_ofreq is not None:
        df_ofreq = df_ofreq.loc[cell_type_lst,:]
    
    df = pd.DataFrame( columns = ['cell_type', 'marker_gene', 'occ_freq_target' ]) 
    cnt = 0
    for ct in cell_type_lst:
        
        b = df_mkr.loc[ct,:] > 0
        genes = df_mkr.columns.values[b]
        
        if df_ofreq is not None:
            # freq = -np.array(df_ofreq.loc[ct,b])

            of1 = df_ofreq.loc[ct,genes].astype(np.float32)
            # of2 = df_ofreq.loc[df_mkr.index.values != ct,b].mean(axis=0).astype(np.float32) + MIN_OFREQ
            of2 = df_ofreq.loc[df_mkr.index.values != ct,genes].max(axis=0).astype(np.float32) + MIN_OFREQ        
            score = np.log2(of1/of2)

            odr = (-score).argsort()
            genes = genes[odr]
            score = score[odr]
        
        
        for k, g in enumerate(genes):
            df.loc[cnt, 'cell_type'] = ct
            df.loc[cnt, 'marker_gene'] = g
            
            if df_ofreq is not None:
                of = df_ofreq[g]
                b = df_mkr.index.values == ct
                of_target = of[ct]
                of_other = of[~b]
                of_other_max = np.max(of_other)
                of_other_mean = np.mean(of_other)

                df.loc[cnt, 'occ_freq_target'] = of_target
                df.loc[cnt, 'occ_freq_non_target'] = of_other_mean
                # df.loc[cnt, 'occ_freq_2nd_max'] = of_other_max

                of_ary = -np.array(of)
                odr = of_ary.argsort()            
                # df.loc[cnt, 'cell_type_2nd_max'] = cell_type_lst[odr[1]]
                # df.loc[cnt, 'score'] = score[k]
                
            cnt += 1
            
    df.to_csv(file)
    return df
 
    
def load_marker_mat( file = 'sc_makers.csv', N_mkrs = 24, score_ref = 'mean' ):
    
    print('Loading markers .. ')
    df = pd.read_csv(file, index_col = 0)
    
    genes_lst = list(df['marker_gene'])
    genes_lst.sort()
    genes_array = np.array(genes_lst)
    
    if score_ref == 'max':
        df['score'] = np.log2(df['occ_freq_target']/(df['occ_freq_2nd_max'] + MIN_OFREQ))
    else:
        df['score'] = np.log2(df['occ_freq_target']/(df['occ_freq_other_mean'] + MIN_OFREQ))
        
    cell_types = list(set(df['cell_type']))
    df_mkr = pd.DataFrame(index = cell_types, columns = genes_array)
    df_mkr.loc[:,:] = 0
    
    for ct in cell_types:
        b = df['cell_type'] == ct
        df_sel = df.loc[b,:]
        scores = np.array(-df_sel['score'])
        odr = scores.argsort()
        if len(odr) <= N_mkrs:
            odr_sel = odr
        else:
            odr_sel = odr[:N_mkrs]
        genes_sel = df_sel.iloc[odr_sel].marker_gene
        df_mkr.loc[ct, genes_sel] = 1
        if len(genes_sel) < N_mkrs:
            print('   %20s has %3i markers only.' % (ct, len(genes_sel)))
    print('Loading done for %i cell types' % len(cell_types))
        
    return df_mkr
        
    
def MkrCnt_Train( X_ref = None, cell_type_ref = None, \
                  mkr_mat = None, N_mkrs = 18, Fold_th = 20, of_th = 0.9, \
                  ttype = 's', stype = 'mean', cell_types_to_excl = [] ):
    
    if mkr_mat is not None:
        df_mkr = mkr_mat.div(mkr_mat.sum(axis=1), axis=0)        
        df_score, cell_type_pred, best_score = get_score( X_ref, df_mkr )
        df_param = get_gmm_param( df_score, cell_type_true = None )        
    else:        
        cell_type_true = np.array(cell_type_ref)
        cell_types_to_test = (list(set(cell_type_ref)))
        if len(cell_types_to_excl) > 0:
            for ct in cell_types_to_excl:
                if ct in cell_types_to_test: 
                    cell_types_to_test.remove(ct)

        ## Generate marker mat
        df_mkr, df_score, cell_type_pred, df_ofreq = \
            gen_marker_profile( X_ref, cell_type_true, N_mkrs, \
                                cell_types_to_excl, of_th = of_th, \
                                ttype = ttype, stype = stype)
        
        df_param = get_gmm_param( df_score, cell_type_true = cell_type_true )

    # b = df_mkr.sum(axis=0) == 1
    # df_mkr = df_mkr.loc[:,b]  
    
    if ttype != 's': df_score = None
    else: 
        df_score['cell_type_ref'] = cell_type_ref
        df_score['cell_type_pred'] = cell_type_pred
        
    return df_mkr, df_param, df_score, df_ofreq


def MkrCnt_Test( X_test, df_mkr_mat, df_param, df_score_ref = None, \
                 min_th = 0.2, log_transformed = True, \
                 ttype = 's', clust_algo = 'gm', N_clusters = None, N_pca = 15, \
                 p_maj = 0.2, p_min = 0.2, cluster_label = None, X_pca = None, \
                 target_FPR = 0.05, verbose = True ):

    df_results = pd.DataFrame( index = X_test.index.values, \
                               columns = ['cluster_label']) 
    
    cells_to_test = df_mkr_mat.index.values
        
    if log_transformed:
        pass
    else:
        X_test = X_test.div((X_test.sum(axis=1)+1e-10)/100000, axis=0)
        X_test = np.log2(1+X_test)
        
    if (cluster_label is None) | (X_pca is None):
            
        if (clust_algo != 'gm') &  (clust_algo != 'km'):
            print('Supported clustering algorithms are gmm and km')
            return
        
        if N_clusters is None:
            N_clusters = int(math.sqrt(X_test.shape[0])/2)
        
        cluster_label, X_pca = \
            pca_and_clustering( X_test, N_clusters, clust_algo = clust_algo, \
                                N_pca_components = N_pca)
        
    N_clusters = len(set(cluster_label))
    if verbose: print_clustering_info(cluster_label)

    df_results['cluster_label'] = cluster_label.astype(str)
    df_score, cell_type_pred, best_score = get_score( X_test, df_mkr_mat ) 

    if isinstance(target_FPR, float):
        target_FPR_lst = [target_FPR]
    else: target_FPR_lst = target_FPR
        
    if df_score_ref is not None:
        cell_type_pred_train = df_score_ref['cell_type_pred']
        cell_type_true_train = df_score_ref['cell_type_ref']
        df_score_train = df_score_ref.drop(columns = ['cell_type_ref', 'cell_type_pred'])
            
    for cnt, tfdr in enumerate(target_FPR_lst):        
        # df_thresh = get_threshold_using_param( df_param, Target_FPR = tfdr )
        if df_score_ref is not None:
            df_thresh = get_threshold( df_score_train, \
                                       cell_type_pred_train, cell_type_true_train, \
                                       Target_FPR = tfdr, verbose = False )                
        else:
            df_thresh = get_threshold_using_param( df_param, Target_FPR = tfdr )
        
        b = df_thresh['threshold'] < min_th
        df_thresh.loc[b,'threshold'] = min_th
        
        b = df_thresh['threshold'] > 0.9
        df_thresh.loc[b,'threshold'] = 0.9
        
        cell_type_pred_c = correct_pred( cell_type_pred, best_score, df_thresh )
        
        cell_type_pred_c = \
            revise_unassinged( cluster_label, cell_type_pred_c, \
                               X_pca = X_pca, df_score = df_score,\
                               p_maj = p_maj, p_min = p_min )

        if len(target_FPR_lst) == 1: cn = 'cell_type_pred'
        else: cn = 'MkrCnt, tFPR_%g' % tfdr            
        df_results[cn] = cell_type_pred_c
            
    return df_results, df_score


def MkrCnt_Ref( X_ref, cell_type_ref, X_test = None, df_mkr_mat = None, \
            N_mkrs = 18, of_th = 0.85, min_th = 0.2, ttype = 's', stype = 'mean', \
            cell_types_to_excl = [], log_transformed = False, \
            clust_algo = 'gm', N_clusters = None, \
            p_maj = 0.2, p_min = 0.2, cluster_label = None, X_pca = None, N_pca = 15, \
            target_FPR = 0.05, file_to_save_marker = None, verbose = True ):

    if ((X_ref is None) | (cell_type_ref is None)):
        print('ERROR:required information  missing .')
        return None
    else:
        if verbose: print('MarkerCount-Ref running ..  ')
    
        X_r = X_ref.copy(deep = True)
        genes_ref = (X_r.columns.values)
        gn_dict = {}
        for gn in list(genes_ref):
            gn_dict[gn] = gn.upper()
        X_r.rename(columns = gn_dict, inplace = True)
        genes_ref = (X_r.columns.values)

        if X_test is not None:
            
            X_t = X_test.copy(deep = True) 
            genes_test = (X_t.columns.values)
            gn_dict = {}
            for gn in list(genes_test):
                gn_dict[gn] = gn.upper()
            X_t.rename(columns = gn_dict, inplace = True)
            genes_test = (X_t.columns.values)

            genes_common = list(set(genes_ref).intersection(genes_test))
            X_r = X_r.loc[:,genes_common]
            X_t = X_t.loc[:,genes_common]
        else:
            genes_common = genes_ref

        if verbose | ( len(genes_test) != len(genes_ref) ):
            print('Ref. has %i genes, Test data has %i genes -> Num. common %i' % (len(genes_ref), len(genes_test), len(genes_common)))
        
        df_mkr, df_param, df_score_ref , df_ofreq = \
            MkrCnt_Train( X_r, cell_type_ref, mkr_mat = df_mkr_mat, \
                          N_mkrs = N_mkrs, of_th = of_th, ttype = ttype, \
                          stype = stype, cell_types_to_excl = cell_types_to_excl )
        if verbose: print_n_mkrs( df_mkr, N_mkrs )
        
        if X_test is not None:
            df_res, df_score = MkrCnt_Test( X_t, df_mkr, df_param, \
                                            df_score_ref = df_score_ref, min_th = min_th, \
                             log_transformed = log_transformed, clust_algo = clust_algo, \
                             N_clusters = N_clusters, N_pca = N_pca, \
                             cluster_label = cluster_label, X_pca = X_pca, \
                             p_maj = p_maj, p_min = p_min, \
                             target_FPR = target_FPR, verbose = verbose )   
        else:
            df_res, df_score = None, None
            
        if (file_to_save_marker is not None) & (isinstance(file_to_save_marker, str)):
            save_markers( df_mkr, file = file_to_save_marker+'_list.csv', df_ofreq = df_ofreq, ct_exc = cell_types_to_excl )
            df_mkr.to_csv(file_to_save_marker+'_matrix.csv')
        
        if verbose: print('MarkerCount-Ref processing done.')
        
    return df_res, df_mkr, df_param, df_score


def MarkerCount_Ref( X_ref, cell_type_ref, X_test = None, df_mkr_mat = None, \
            N_mkrs = 18, of_th = 0.85, cell_types_to_excl = [], log_transformed = False, \
            N_clusters = None, cluster_label = None, X_pca = None, N_pca = 15, \
            target_FPR = 0.05, file_to_save_marker = None, verbose = True ):

    df_res, df_mkr, df_param, df_score = \
        MkrCnt_Ref( X_ref = X_ref, cell_type_ref = cell_type_ref, X_test = X_test, \
            df_mkr_mat = df_mkr_mat, N_mkrs = N_mkrs, of_th = of_th, \
            cell_types_to_excl = cell_types_to_excl, log_transformed = log_transformed, \
            N_clusters = N_clusters, cluster_label = cluster_label, X_pca = X_pca, N_pca = N_pca, \
            target_FPR = target_FPR, file_to_save_marker = file_to_save_marker, verbose = verbose )

    if X_test is not None:
        return df_res
    else:
        return df_mkr


def get_best_n_markers( df_ofreq, ct, of_th = 0.9, N_mkr = 0, stype = 'mean', mkrs_known = None):
    
    genes = df_ofreq.columns.values

    b = df_ofreq.index.values == ct
    of1 = df_ofreq.loc[ct,:].astype(np.float32) + MIN_OFREQ
    if stype == 'max':
        of2 = df_ofreq.loc[~b,:].max(axis=0).astype(np.float32) + MIN_OFREQ
    else:
        of2 = df_ofreq.loc[~b,:].mean(axis=0) + MIN_OFREQ

    score = (of1)*((1-of2)**2)

    if N_mkr == 0:
        bs1 = (of1 >= of_th) & (of2 < (1-of_th))
        return genes[bs1]
    else:
        bf1 = (of1 >= of_th) & (of2 < (1-of_th))  
        s = np.array(score)
        b = np.full(len(s), False)
        odr = (-s).argsort()
        b[odr[:N_mkr]] = True
        if mkrs_known is None:
            bs1 = b | bf1
            return genes[bs1]
        else:
            for k in range(1,N_mkr):
                    gn = list( set(mkrs_known).union(genes[odr[:k]]) )
                    if len(gn) >= N_mkr:
                        break
            return gn


def MkrCnt( X, mkr_mat, min_th = 0.2, N_mkrs_max = 18, N_mkrs_min = 3, of_th = 0.85, \
            target_FPR = 0.12, clustering = 'gm', N_clusters = None, \
            cluster_label = None, X_pca = None, N_pca = 15, \
            p_maj = 0.2, p_min = 0.2, init_fpr = 0.08, \
            log_transformed = False, Loop = 1, sr = 0.8, \
            ct_map = None, verbose = False ):

    N_mkrs = N_mkrs_max
    N_mkrs_min = min(N_mkrs_min, N_mkrs_max)
    
    df_comp = pd.DataFrame( index = X.index.values ) 
    cells_to_test = list(mkr_mat.index.values)

    if verbose: print('MarkerCount running ..  ')
    
    ######## Preprocessing for test data
    if log_transformed:
        X_for_clustering = X
    else:
        X = X.div((X.sum(axis=1) + 1e-10) / 100000, axis=0)
        X = np.log2(1 + X)
        X_for_clustering = X

    if (cluster_label is None) | (X_pca is None):
        if N_clusters is None:
            N_clusters = int(math.sqrt(X.shape[0])/2)    
        cluster_label_tmp, X_pca = pca_and_clustering(X_for_clustering, N_clusters, clust_algo = clustering, \
                                   N_pca_components = N_pca)
        if (cluster_label is None):
            cluster_label = cluster_label_tmp
    else:
        N_clusters = len(set(cluster_label))
        
    df_comp['cluster_label'] = cluster_label
    if verbose: print_clustering_info(cluster_label)

    df_mkr = mkr_mat
    df_score, cell_type_pred, best_score = get_score( X, df_mkr )
    df_param = get_gmm_param( df_score, cell_type_true = None )
    
    if (Loop > 0): 
        
        df_thresh = get_threshold_using_param( df_param, Target_FPR = init_fpr )

        b = df_thresh['threshold'] < min_th
        df_thresh.loc[b, 'threshold'] = min_th

        cell_type_pred = correct_pred( cell_type_pred, best_score, df_thresh )
    
        df_mkr2 = df_mkr.copy(deep = True)
        ct_pred2 = copy.deepcopy(cell_type_pred)
        cluster_lst = list(set(cluster_label))
        cluster_lst.sort()
        df_score2 = df_score.copy(deep = True)
        best_score2 = copy.deepcopy(best_score)
        
        odr = np.array(best_score2).argsort()
        idx = int(len(odr)*sr)
        sth = best_score2[odr[idx]]

        for k in range(Loop):

            for n, clst in enumerate(cluster_lst):
                b = cluster_label == clst
                ct_sel = ct_pred2[b]
                ct_lst_sorted, cell_num = get_type_num(ct_sel, verbose = False)
                total_num = np.sum(cell_num)
                for m, ct in enumerate(ct_lst_sorted):
                    if (ct != STR_UNASSIGNED) :
                        b1 = (ct_pred2 == ct)
                        ## if the cell type comprises at least 10% of a cluster
                        if (cell_num[m] >= total_num*0.1) & (cell_num[m] >= 10): 
                            pass
                        else:
                            ct_pred2[b&b1] = STR_UNASSIGNED

            df_ofreq_tmp_all = get_occurrence_freq(X, ct_pred2, cells_to_test)

            for ct in cells_to_test:
                b = df_mkr.loc[ct,:] > 0
                ct_mkrs = df_mkr.columns.values[b]
                ct_mkrs = set(ct_mkrs).intersection(X.columns.values)
                
                if len(ct_mkrs) > N_mkrs:
                     if ct in list(df_ofreq_tmp_all.index.values):
                        mkrs = get_best_n_markers( df_ofreq_tmp_all.loc[:,ct_mkrs], ct, of_th = of_th, \
                                                   N_mkr = N_mkrs)
                        df_mkr2.loc[ct,:] = 0
                        df_mkr2.loc[ct,mkrs] = 1
                        
                elif len(ct_mkrs) < N_mkrs_min:                    
                    if ct in list(df_ofreq_tmp_all.index.values):
                        mkrs = get_best_n_markers( df_ofreq_tmp_all, ct, of_th = of_th, \
                                                   N_mkr = N_mkrs_min, mkrs_known = ct_mkrs)
                        df_mkr2.loc[ct,:] = 0
                        for mkr in mkrs:
                            if mkr in list(df_mkr2.columns.values):
                                df_mkr2.loc[ct,mkr] = 1
                            else:
                                df_mkr2[mkr] = 0
                                df_mkr2.loc[ct, mkr] = 1

            df_score2, ct_pred2, best_score2 = get_score( X, df_mkr2 )
            ct_pred = ct_pred2

            if verbose: 
                print_n_mkrs( df_mkr2, N_mkrs )
                if ct_map is not None:
                    ct_pred_tmp = rename_cell_type(ct_pred, ct_map)
                    print( 'Loop %i: %5i/%5i' % \
                       (k+1, np.sum( ct_pred_tmp == cell_type), len(cell_type)))

        cell_type_pred = ct_pred
        df_score = df_score2
        df_param = get_gmm_param( df_score, cell_type_true = cell_type_pred )

    ## Loop end
    
    if isinstance(target_FPR, float):
        target_FPR_lst = [target_FPR]
    else:
        target_FPR_lst = target_FPR
        
    df_thresh_all = pd.DataFrame(index = df_thresh.index.values)
    for cnt, tfdr in enumerate(target_FPR_lst):
        df_thresh = get_threshold_using_param( df_param, Target_FPR = tfdr )
        
        b = df_thresh['threshold'] < min_th
        df_thresh.loc[b, 'threshold'] = min_th
    
        b = df_thresh['threshold'] > 0.9
        df_thresh.loc[b, 'threshold'] = 0.9
    
        cell_type_pred_c = correct_pred( cell_type_pred, best_score, df_thresh )
        if cluster_label is not None:
            cell_type_pred_c = \
                revise_unassinged( cluster_label, cell_type_pred_c, X_pca = X_pca, \
                                   df_score = df_score, p_maj = p_maj, p_min = p_min )
    
        if ct_map is not None:
            cell_type_pred_c = rename_cell_type(cell_type_pred_c, ct_map)
            
        if len(target_FPR_lst) == 1:
            cn = 'cell_type_pred'
        else:
            cn = 'Target_FPR_%g' % tfdr
        df_comp[cn] = cell_type_pred_c
        df_thresh_all[cn] = df_thresh['threshold']
        
        if verbose: print('MarkerCount processing done.')
        
    return df_comp, df_thresh_all, df_score, df_param


def MarkerCount( X, mkr_mat, target_FPR = 0.12, \
            N_clusters = None, cluster_label = None, X_pca = None, N_pca = 15, \
            log_transformed = False, verbose = False ):

    df_pred, df_thresh_all, df_score, df_param = \
        MkrCnt( X, mkr_mat, target_FPR = target_FPR, \
            N_clusters = N_clusters, cluster_label = cluster_label, X_pca = X_pca, N_pca = N_pca, \
            log_transformed = log_transformed, verbose = verbose )

    return df_pred
