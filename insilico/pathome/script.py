#!/usr/bin/env python
# coding: utf-8




import itertools
import gzip
import os
import glob
import json
import scipy
import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from sklearn.feature_extraction.text import CountVectorizer
import os, zipfile, csv, codecs, gzip, itertools, json, math, threading, requests, scipy
from datetime import datetime
from pprint import pprint
from sklearn.decomposition import LatentDirichletAllocation





def load_lrp_document_term(tissue, matrix):
    return scipy.sparse.load_npz('lrp_document_term.'+tissue.lower()+'.'+matrix+'.npz'); print(lrp_document_term.shape)

def load_lrp_document_info(tissue):
    return pd.read_csv('lrp_document_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_document_info.shape)

def load_lrp_term_info(tissue):
    return pd.read_csv('lrp_term_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_term_info.shape)

def load_lrp_lda_components(tissue):
    return scipy.sparse.load_npz('lrp_document_term.'+tissue.lower()+'.lda_components.npz'); print(lrp_document_term_lda_components.shape)





def load_mean_age(tissue):
    return pd.read_csv('lrp_term_info.'+tissue+'.mean_age.tsv', sep='\t')['mean_age']

def load_cox_coef(tissue):
    return pd.read_csv('lrp_term_info.'+tissue+'.cox_coef.tsv', sep='\t')['cox_coef']

def load_lda_topic(tissue):
    return pd.read_csv('lrp_term_info.'+tissue+'.lda_topic.tsv', sep='\t')['lda_topic']

def load_term(tissue):
    return pd.read_csv('lrp_term_info.'+tissue+'.term.tsv', sep='\t')





def run_age_aggregated_tfidf():
    lrp_document_term = load_lrp_document_term(tissue='t28', matrix='onehot')
    lrp_document_term_tfidf = TfidfTransformer(smooth_idf=False, norm=None).fit_transform(lrp_document_term)
    lrp_document_info = load_lrp_document_info(tissue='t28')    
    lrp_document_info['T3_years'] = (lrp_document_info.T3/365.25).round().astype(int)
    surv_ages = lrp_document_info.T3_years.tolist()
    ages, counts = np.unique(lrp_document_info.T3_years, return_counts=True)    
    X_tfidf_filter_compressed = np.zeros((max(ages)+1,lrp_document_term_tfidf.shape[1])); print(X_tfidf_filter_compressed.shape)
    X_binary_filter_compressed = np.zeros((max(ages)+1,lrp_document_term.shape[1])); print(X_binary_filter_compressed.shape)
    for i, v in enumerate(surv_ages):
        X_tfidf_filter_compressed[v,:]=X_tfidf_filter_compressed[v,:]+lrp_document_term_tfidf[i,:].toarray()
        X_binary_filter_compressed[v,:]=X_binary_filter_compressed[v,:]+lrp_document_term[i,:].toarray()
    for i, v in enumerate(ages):
        X_tfidf_filter_compressed[v,:]=X_tfidf_filter_compressed[v,:]/counts[i]
        X_binary_filter_compressed[v,:]=X_binary_filter_compressed[v,:]/counts[i]  
    scipy.sparse.save_npz('lrp_document_term.'+tissue.lower()+'.compressed_tfidf.npz', scipy.sparse.csr_matrix(X_tfidf_filter_compressed))





def run_lda(tissue):
    lrp_document_term = load_lrp_document_term(tissue='t28', matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue='t28')
    for i in range(2, 150):
        lda = LatentDirichletAllocation(n_components=i)
        lda.fit(lrp_document_term)
        X_fit_transform = lda.fit_transform(lrp_document_term)
        dump(lda, tissue + '_lda_' + str(i) + '.joblib')
        dump(X_fit_transform, tissue + '_lda_' + str(i) + '_X_fit_transform.joblib')





def run_distance_matrix(tissue, matrix, metric):
    lrp_document_term = load_lrp_document_term(tissue, matrix)
    dist = pairwise_distances(lrp_document_term.transpose(), metric=metric)    
    np.save('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.npy', dist)    
    
def load_distance_matrix(tissue, matrix, metric):
    return np.load('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.npy')    





def run_umap_embedding(tissue, matrix, metric):
    dist = load_distance_matrix(tissue, matrix, metric) 
    umap_embedding = umap.UMAP(random_state=0, metric='precomputed').fit_transform(dist)
    np.save('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.umap.npy', umap_embedding)            
    
def load_umap_embedding(tissue, matrix, metric):
    return np.load('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.umap.npy')            





def run_tsne_embedding(tissue, matrix, metric):
    dist = load_distance_matrix(tissue, matrix, metric) 
    tsne_embedding = TSNE(n_components=2, random_state=0, metric='precomputed').fit_transform(dist)        
    np.save('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.tsne.npy', tsne_embedding)
    
def load_tsne_embedding(tissue, matrix, metric):
    return np.load('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.umap.npy')





def load_embedding(tissue, matrix, metric, method):
    return np.load('lrp_document_term.'+tissue.lower()+'.'+matrix+'.dist.'+metric+'.'+method+'.npy')            





def run_pca(tissue, matrix):
    lrp_document_term = load_lrp_document_term(tissue, matrix)
    pca = PCA(n_components=2, random_state=0)
    pca_decomposition = pca.fit_transform(lrp_document_term.todense())
    np.save('lrp_document_term.'+tissue+'.'+matrix+'.pca.npy', pca_decomposition)
    dump(pca, 'lrp_document_term.'+tissue+'.'+matrix+'.pca.joblib')
    
def load_pca(tissue, matrix):
    return np.load('lrp_document_term.'+tissue+'.'+matrix+'.pca.npy')





def plot_pca(tissue, age_start, age_end):
    my_font_size=30
    matrix='compressed_tfidf'
    embedding = load_pca(tissue=tissue, matrix=matrix)
    c = load_mean_age(tissue)
    pca = load('lrp_document_term.'+tissue+'.'+matrix+'.pca.joblib')
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_title('PCA', fontsize=my_font_size)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    sc = ax.scatter(embedding[:,0][age_start:age_end], embedding[:,1][age_start:age_end], cmap='coolwarm', c=np.arange(embedding.shape[0])[age_start:age_end], s=500)
    ax.set_xlabel('PC1 - ' + str(round(100*pca.explained_variance_ratio_[0],2)) + '%', fontsize=my_font_size)
    ax.set_ylabel('PC2 - ' + str(round(100*pca.explained_variance_ratio_[1],2)) + '%', fontsize=my_font_size)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines['right'].set_visible(False); ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False); ax.spines['top'].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.ax.tick_params(labelsize=my_font_size)
    cbar.set_label('Age', fontsize=my_font_size)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(size=0)
    plt.savefig(tissue+'.'+'compressed_tfidf'+'.'+'pca'+'.age.'+str(age_start)+'.'+str(age_end)+'.png', dpi=1200, transparent=True)
    plt.close(fig)

def plot_embedding_mean_age(tissue, matrix, metric, method):
    my_font_size=30    
    embedding = load_embedding(tissue=tissue, matrix=matrix, metric=metric, method=method)
    c = load_mean_age(tissue)
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_title('t-SNE', fontsize=my_font_size) 
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    sc = plt.scatter(embedding[:,0], embedding[:,1], c=c, cmap='coolwarm', norm=colors.TwoSlopeNorm(vmin=c.min(), vcenter=c.mean(), vmax=c.max()), s=1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines['right'].set_visible(False); ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False); ax.spines['top'].set_visible(False)    
    cbar = fig.colorbar(sc, ax=ax)
    cbar.ax.tick_params(labelsize=my_font_size)    
    cbar.set_label('Mean age', fontsize=my_font_size)    
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(size=0)    
    plt.savefig(tissue+'.'+matrix+'.'+metric+'.'+method+'.mean_age.png', dpi=1200, transparent=True)
    plt.close(fig)
    
def plot_embedding_cox_coef(tissue, matrix, metric, method):
    embedding = load_embedding(tissue=tissue, matrix=matrix, metric=metric, method=method)
    c = load_cox_coef(tissue)    
    fig, ax = plt.subplots(figsize=(5,5))
    ax.axes.xaxis.set_ticks([]); ax.axes.yaxis.set_ticks([]); ax.axis('off')
    sc = plt.scatter(embedding[:,0], embedding[:,1], c=c, cmap='coolwarm', norm=colors.TwoSlopeNorm(vmin=c.min(), vcenter=0, vmax=c.max()), s=1)
    cbar = fig.colorbar(sc, ax=ax, extend='both')
    cbar.outline.set_visible(False)
    plt.savefig(tissue+'.'+matrix+'.'+metric+'.'+method+'.cox_coef.png', dpi=300, transparent=True)
    plt.close(fig)
    
def plot_embedding_lda_topic(tissue, matrix, metric, method, topics):
    embedding = load_embedding(tissue=tissue, matrix=matrix, metric=metric, method=method)
    c = load_lda_topic(tissue)    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.axes.xaxis.set_ticks([]); ax.axes.yaxis.set_ticks([]); ax.axis('off')
    sc = plt.scatter(embedding[:,0], embedding[:,1], c=c, cmap='tab20', s=10)
    plt.savefig(tissue+'.'+matrix+'.'+metric+'.'+method+'.lda_topic.png', dpi=300, transparent=True)
    plt.close(fig)
    if topics==True:
        for topic in range(c.max()+1):
            fig, ax = plt.subplots(figsize=(10,10))
            ax.axes.xaxis.set_ticks([]); ax.axes.yaxis.set_ticks([]); ax.axis('off')
            sc = plt.scatter(embedding[:,0], embedding[:,1], c=(c==topic), cmap='coolwarm', s=10)
            plt.title('Topic '+str(topic))
            plt.savefig(tissue+'.'+matrix+'.'+metric+'.'+method+'.lda_topic.topic'+str(topic)+'.png', dpi=300, transparent=True)
            plt.close(fig)
            plt.show()





def load_regression(tissue, model_type):
    return load('lrp_document_term.'+tissue+'.'+model_type+'.joblib')

def save_regression(tissue, model_type, model):
    dump(model, 'lrp_document_term.'+tissue+'.'+model_type+'.joblib')

def run_linear_regression(tissue, matrix):
    lrp_document_term = load_lrp_document_term(tissue=tissue, matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    model = LinearRegression().fit(X_train, y_train)
    save_regression(tissue=tissue, model_type='LinearRegression', model=model)
    
def run_mpl_regression(tissue, matrix, max_iter):
    lrp_document_term = load_lrp_document_term(tissue=tissue, matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    model = MLPRegressor(random_state=0, max_iter=max_iter).fit(X_train, y_train)
    save_regression(tissue=tissue, model_type='MLPRegressor'+str(max_iter), model=model)
    
def plot_regression_per_topic(tissue, model_type):
    model = load_regression(tissue=tissue, model_type=model_type)
    lrp_document_term = load_lrp_document_term(tissue=tissue, matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    tt = lrp_document_info['topic'].tolist()
    xx = y/365.25
    yy = model.predict(X)/365.25
    df = pd.DataFrame(zip(xx, yy, tt), columns=['x', 'y', 'topic'])
    sns.lmplot(x='x', y='y', hue='topic', data=df, scatter=False, legend=False, ci=False, line_kws={'color':'black', 'linewidth':1})
    plt.savefig(tissue+'.'+model_type+'.age_prediction_per_topic.png', dpi=300)
    plt.title('Age prediction per topic')    
    plt.close()
    
def plot_regression(tissue, model_type):
    model = load_regression(tissue=tissue, model_type=model_type)
    lrp_document_term = load_lrp_document_term(tissue=tissue, matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    xx = np.digitize(y_test/365.25, bins=np.arange(0,100,10))*10-5
    yy = model.predict(X_test)/365.25
    df = pd.DataFrame(zip(xx, yy), columns=['Age', 'Predicted age'])
    sns.boxplot(x='Age', y='Predicted age', data=df, showfliers=False, color='grey')  
    plt.title('Age prediction (R²='+str(np.round(model.score(X_test, y_test), 2))+')')
    plt.savefig(tissue+'.'+model_type+'.age_prediction.png', dpi=600, transparent=True)
    plt.close()





def run_permutation_importance(tissue, model_type, max_samples):
    model = load('lrp_document_term.'+tissue+'.'+model_type+'.joblib')

    lrp_document_term = scipy.sparse.load_npz('lrp_document_term.'+tissue.lower()+'.onehot.npz'); print(lrp_document_term.shape)
    lrp_document_info = pd.read_csv('lrp_document_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_document_info.shape)
    lrp_term_info = pd.read_csv('lrp_term_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_term_info.shape)

    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    
    result = permutation_importance(model, X_train.todense(), y_train, n_jobs=8, n_repeats=5, max_samples=max_samples)
    dump(result, 'lrp_document_term.'+tissue+'.'+model_type+'.pfi.result')
    
def plot_permutation_importance(tissue, model_type):
    result =  load('lrp_document_term.'+tissue+'.'+model_type+'.pfi.result')
    lrp_term_info = pd.read_csv('lrp_term_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_term_info.shape)
    
    feature_names = lrp_term_info['feature']
    perm_sorted_idx = np.flip(result.importances_mean.argsort())
    fig, ax = plt.subplots()
    ax.boxplot(
        result.importances[perm_sorted_idx][0:20].T,
        vert=False,
        labels=feature_names[perm_sorted_idx][0:20],
    )
    plt.title('Permutation Feature Importance')
    plt.xlabel(r'$\Delta$'+'R²')

    plt.tight_layout()
    plt.savefig(tissue+'.'+model_type+'.plot_permutation_importance.png', dpi=600, transparent=True)    
    plt.close()





def run_permutation_importance_per_topic(tissue, model_type, n_repeats, n_clust, n_samples):
    model = load('lrp_document_term.'+tissue+'.'+model_type+'.joblib')

    lrp_document_term = scipy.sparse.load_npz('lrp_document_term.'+tissue.lower()+'.onehot.npz'); print(lrp_document_term.shape)
    lrp_document_info = pd.read_csv('lrp_document_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_document_info.shape)
    lrp_term_info = pd.read_csv('lrp_term_info.'+tissue.lower()+'.tsv', sep='\t'); print(lrp_term_info.shape)

    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    
    clust_labels = lrp_term_info['topic'].tolist(); print(len(clust_labels))
    X_perm = X_train.copy(); print(X_perm.shape)
    y_perm = y_train.copy(); print(y_perm.shape)

    cluster_scores = np.zeros((n_clust, n_repeats))
    importances = np.zeros((n_clust, n_repeats))
    importances_mean = np.zeros((n_clust, 1))
    importances_std = np.zeros((n_clust, 1))
    X_perm_lil_orig = X_perm.tolil()
    X_perm_csc_orig = X_perm.tocsc()

    for repeat_num in range(n_repeats):
        print(repeat_num)
        for clust_num in range(n_clust):
            clust_idxs = (np.array(clust_labels)==clust_num).nonzero()[0]
            ###clust_idxs = feature_info[feature_info['topic']==clust_num].sort_values(by='topic_max', ascending=False)[0:10].index.to_numpy().tolist()
            samples = np.random.randint(0, X_perm.shape[0], n_samples)        
            X_perm_lil = X_perm_lil_orig[samples,:].copy()
            X_perm_csc = X_perm_csc_orig[samples,:].copy()
            for clust_idx in clust_idxs:
                index = np.arange(X_perm_lil.shape[0]); np.random.shuffle(index)
                X_perm_lil[:, clust_idx] = X_perm_csc[index, clust_idx]
            X_perm_csr = X_perm_lil.tocsr()
            cluster_scores[clust_num, repeat_num] = model.score(X_perm_csr, y_perm[samples])
    importances = model.score(X_perm, y_perm) - cluster_scores 
    importances_mean = importances.mean(axis=1)
    importances_std = importances.std(axis=1)

    np.save('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances', importances)
    np.save('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances_mean', importances_mean)
    np.save('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances_std', importances_std)
    
def plot_permutation_importance_per_topic(tissue, model_type, topic):
    importances = np.load('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances.npy')
    importances_mean = np.load('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances_mean.npy')
    importances_std = np.load('lrp_document_term.'+tissue+'.'+model_type+'.topic_pfi.importances_std.npy')    
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    feature_names = np.arange(lrp_document_info['topic'].max()+1)
    perm_sorted_idx = importances_mean.argsort()
    fig, ax = plt.subplots(figsize=(10, 10))
    bplot = ax.boxplot(
        importances[perm_sorted_idx].T,
        vert=False,
        patch_artist=True,
        labels=['Topic '+str(i) for i in feature_names[perm_sorted_idx]],
    )
    bplot['boxes'][feature_names[perm_sorted_idx].tolist().index(topic)].set_facecolor('red')
    plt.title('Permutation Topic Importance')
    plt.xlabel(r'$\Delta$'+'R²')
    plt.tight_layout()

    plt.savefig(tissue+'.'+model_type+'.plot_permutation_importance_per_topic.png', dpi=600, transparent=True)
    plt.close()    





def run_topic_slopes(tissue, model_type):
    model = load('lrp_document_term.'+tissue+'.'+model_type+'.joblib')
    lrp_document_term = load_lrp_document_term(tissue=tissue, matrix='onehot')
    lrp_document_info = load_lrp_document_info(tissue=tissue)
    X = lrp_document_term.copy()
    y = lrp_document_info['T3'].to_numpy().copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)    
    n_clusts = lrp_document_info['topic'].max()+1
    clust_num = 0
    cluster_slopes = np.zeros((n_clusts, 1))
    clust_labels = lrp_document_info['topic'].tolist()
    for clust_num in range(n_clusts):
        clust_idxs = (np.array(clust_labels)==clust_num).nonzero()[0]
        x_slope = y[clust_idxs]
        y_slope = model.predict(X)[clust_idxs]
        m, b = np.polyfit(x_slope, y_slope, 1)
        cluster_slopes[clust_num, 0] = m
    dump(cluster_slopes, 'lrp_document_term.'+tissue+'.'+model_type+'.topic_slopes.result')
        
def plot_topic_slopes(tissue, model_type, topic):
    df =  pd.DataFrame(load('lrp_document_term.'+tissue+'.'+model_type+'.topic_slopes.result'), columns=['slope'])
    df = df.sort_values(by='slope')
    df['topic']=df.index
    df['color']=(df['topic']==topic)    
    df['topic']='Topic ' + df['topic'].astype(str)
    fig, ax = plt.subplots(figsize=(20,10))
    sns.pointplot(x='topic', y='slope', hue='color', data=df, join=False)
    plt.legend([],[], frameon=False)
    ax.tick_params(axis='x', rotation=90)
    plt.tight_layout()
    plt.title('Topic age regression slopes')
    plt.savefig(tissue+'.'+model_type+'.topic_slopes.png', dpi=600, transparent=True)
    plt.close()    





def plot_cox(topic):
    df = load_cox_coef(tissue='t28')
    df = df.iloc[2:].reset_index(drop=True)
    df = df.sort_values(by='coef')
    df['topic']=df.index
    df['color']=(df['topic']==topic)        
    df['topic']='Topic' + df['topic'].astype(str)
    plt.subplots(figsize=(10,15))
    sns.pointplot(x='coef', y='topic', hue='color', data=df, join=False)
    plt.legend([],[], frameon=False)
    plt.savefig(tissue+'.topic_cox_coef.png', dpi=300)
    plt.close()    





def init_pubmed():
    xml_files = glob.glob('pubmed/baseline/pubmed21n*.xml.gz'); len(xml_files)
    for xml_file in xml_files:
        print(xml_file)
        xml_file_base = os.path.basename(xml_file).split('.')[0]
        with gzip.open(xml_file, 'r') as infile:
            tree = ET.parse(infile)
            root = tree.getroot()
            PubmedArticles = root.findall('PubmedArticle')
            docs = []
            for PubmedArticle in PubmedArticles:
                doc = {'PMID':None,'Title':None,'Abstract':None,'Year':None,'Chemicals':[],'Mesh':[]}
                MedlineCitation = PubmedArticle.find('MedlineCitation')
                if MedlineCitation is not None:
                    PMID = MedlineCitation.find('PMID')
                    if PMID is not None:
                        doc['PMID'] = PMID.text
                    DateCompleted = MedlineCitation.find('DateCompleted')
                    if DateCompleted is not None:
                        Year = DateCompleted.find('Year')
                        if Year is not None:
                            doc['Year'] = Year.text
                    Article = MedlineCitation.find('Article')
                    if Article is not None:
                        ArticleTitle = Article.find('ArticleTitle')
                        if ArticleTitle is not None:
                            doc['Title'] = ArticleTitle.text
                        Abstract = Article.find('Abstract')
                        if Abstract is not None:
                            AbstractText = Abstract.find('AbstractText')
                            if AbstractText is not None:
                                doc['Abstract'] = AbstractText.text
                    ChemicalList = MedlineCitation.find('ChemicalList')
                    if ChemicalList is not None:
                        Chemicals = ChemicalList.findall('Chemical')
                        if Chemicals is not None:
                            for Chemical in Chemicals:
                                NameOfSubstance = Chemical.find('NameOfSubstance')
                                if NameOfSubstance is not None:
                                    doc['Chemicals'].append(NameOfSubstance.get('UI'))
                    MeshHeadingList = MedlineCitation.find('MeshHeadingList')
                    if MeshHeadingList is not None:            
                        MeshHeadings = MeshHeadingList.findall('MeshHeading')
                        if MeshHeadings is not None:
                            for MeshHeading in MeshHeadings:
                                MeshHeadingChildren = list(MeshHeading)
                                if MeshHeadingChildren is not None:
                                    for MeshHeadingChild in MeshHeadingChildren:
                                        doc['Mesh'].append(MeshHeadingChild.get('UI'))
                    docs.append(doc)
        with gzip.open('xml_file_base + '.json.gz', 'wt') as outfile:
            for doc in docs:
                outfile.write(json.dumps(doc) + '\n')            

    files = glob.glob('pubmed21n*.json.gz'); print(len(files))
    with gzip.open('pubmed21.json.gz', 'wt') as outfile:
        for file in files:
            with gzip.open(file, 'rt') as infile:
                for line in infile:
                    doc = json.loads(line)
                    if doc['PMID'] is None:
                        print(doc)
                    if doc['Year'] is None:
                        doc['Year'] = ''
                    if doc['Abstract'] is None:
                        doc['Abstract'] = ''
                    if doc['Title'] is None:
                        doc['Title'] = ''
                    doc['TitleAndAbstract'] = doc['Title'] + ' ' + doc['Abstract']
                    doc['TitleAndAbstract'] = doc['TitleAndAbstract'].strip()
                    outfile.write(json.dumps(doc) + '\n')

    sct_en_term_info = pd.read_pickle('data/sct_en_term_info.pkl'); print(sct_en_term_info.shape)
    sct_en_term_info
                       
    def read_pubmed():
        with gzip.open('pubmed21.json.gz', 'rt') as infile:
            for line in infile:
                yield (json.loads(line)['TitleAndAbstract'])

    pubmed_term_info = sct_en_term_info.copy()
    pubmed_term_info = pubmed_term_info.drop_duplicates(subset='term')
    pubmed_term_info = pubmed_term_info.reset_index(drop=True)
    pubmed_term_info

    count_vectorizer = CountVectorizer(lowercase=True, binary=True, dtype=np.int8, ngram_range=(1, 4), vocabulary=pubmed_term_info['term'].tolist())
    pubmed_document_term = count_vectorizer.fit_transform(read_pubmed()); print(pubmed_document_term.shape)
    scipy.sparse.save_npz('data/pubmed_sct_en_document_term.npz', pubmed_document_term[:,pubmed_document_term.sum(axis=0).A1>0], compressed=True)

    pubmed_term_info['pubmed_df'] = pubmed_document_term.sum(axis=0).A1
    pubmed_term_info = pubmed_term_info[pubmed_term_info['pubmed_df']>0]
    pubmed_term_info = pubmed_term_info.reset_index(drop=True)
    pubmed_term_info['pubmed_log10df'] = np.log10(pubmed_term_info['pubmed_df'])
    pubmed_term_info['pubmed_idx'] = pubmed_term_info.index
    pubmed_term_info.to_pickle('data/pubmed_sct_en_term_info.pkl')





def init_pubchem():
    CID_PMID = pd.read_csv('pubchem/Compound/Extras/CID-PMID.gz', delimiter='\t', header=None); print(CID_PMID.shape)
    CID_PMID.columns = ['CID', 'PMID', 'type']
    CID_PMID = CID_PMID.drop(['type'], axis=1)
    CID_PMID = CID_PMID.drop_duplicates(subset=['CID', 'PMID'])
    CID_PMID = CID_PMID.reset_index(drop=True)
    CID_PMID.to_pickle('CID_PMID.pkl')    
    
    CID_PMID = pd.read_pickle('CID_PMID.pkl')
    CID_CIDidx = CID_PMID['CID']
    CID_CIDidx = CID_CIDidx.drop_duplicates()
    CID_CIDidx = CID_CIDidx.sort_values()
    CID_CIDidx = CID_CIDidx.reset_index(drop=True)
    CID_CIDidx = CID_CIDidx.to_frame()
    CID_CIDidx['CIDidx'] = CID_CIDidx.index
    CID_CIDidx.to_pickle('CID_CIDidx.pkl')    
    
    CID_PMID = pd.read_pickle('CID_PMID.pkl')
    PMID_PMIDidx = CID_PMID['PMID']
    PMID_PMIDidx = PMID_PMIDidx.drop_duplicates()
    PMID_PMIDidx = PMID_PMIDidx.sort_values()
    PMID_PMIDidx = PMID_PMIDidx.reset_index(drop=True)
    PMID_PMIDidx = PMID_PMIDidx.to_frame()
    PMID_PMIDidx['PMIDidx'] = PMID_PMIDidx.index
    PMID_PMIDidx.to_pickle('PMID_PMIDidx.pkl')    

    CID_PMID = pd.read_pickle('CID_PMID.pkl'); print(CID_PMID.shape)
    CID_CIDidx = pd.read_pickle('CID_CIDidx.pkl'); print(CID_CIDidx.shape)
    PMID_PMIDidx = pd.read_pickle('PMID_PMIDidx.pkl'); print(PMID_PMIDidx.shape)    
    
    CID_PMID_CIDidx = CID_PMID.merge(CID_CIDidx, on='CID', how='left')
    CID_PMID_CIDidx.to_pickle('CID_PMID_CIDidx.pkl')    
    
    CID_PMID_PMIDidx = CID_PMID.merge(PMID_PMIDidx, on='PMID', how='left')
    CID_PMID_PMIDidx.to_pickle('CID_PMID_PMIDidx.pkl')    
    
    CID_CIDidx = pd.read_pickle('CID_CIDidx.pkl')
    PMID_PMIDidx = pd.read_pickle('PMID_PMIDidx.pkl')
    CID_PMID_CIDidx = pd.read_pickle('CID_PMID_CIDidx.pkl')
    CID_PMID_PMIDidx = pd.read_pickle('CID_PMID_PMIDidx.pkl')    
    
    CID_PMID_coo = coo_matrix(([1]*CID_PMID_CIDidx.shape[0], (CID_PMID_PMIDidx['PMIDidx'], CID_PMID_CIDidx['CIDidx'])), 
                              shape=(PMID_PMIDidx.shape[0], CID_CIDidx.shape[0]), dtype=np.int8)
    scipy.sparse.save_npz('CID_PMID_coo.npz', CID_PMID_coo, compressed=True)
    pubchem_document_term = CID_PMID_coo.tocsr()
    pubchem_document_term[pubchem_document_term!=0]=1
    scipy.sparse.save_npz('data/pubchem_document_term.npz', pubchem_document_term, compressed=True)    
    
    CID_MeSH = pd.read_csv('pubchem/Compound/Extras/CID-MeSH', sep='^([0-9]*)\t', header=None, engine='python').loc[:,1:2]; print(CID_MeSH.shape)
    CID_MeSH.columns = ['CID', 'MeSH']
    CID_MeSH.to_pickle('CID_MeSH.pkl')    
    
    CID_CIDidx = pd.read_pickle('CID_CIDidx.pkl')
    CID_MeSH = pd.read_pickle('CID_MeSH.pkl')
    CID_CIDidx_MeSH = CID_CIDidx.merge(CID_MeSH, on='CID', how='left')
    CID_CIDidx_MeSH.to_pickle('CID_CIDidx_MeSH.pkl')    
    
    CID_Title = pd.read_csv('pubchem/Compound/Extras/CID-Title.gz', sep='\t', header=None, encoding='iso-8859-1'); print(CID_Title.shape)
    CID_Title.columns = ['CID', 'Title']
    CID_Title.to_pickle('CID_Title.pkl')    
    
    CID_CIDidx = pd.read_pickle('CID_CIDidx.pkl')
    CID_Title = pd.read_pickle('CID_Title.pkl')
    CID_CIDidx_Title = CID_CIDidx.merge(CID_Title, on='CID', how='left')
    CID_CIDidx_Title.to_pickle('CID_CIDidx_Title.pkl')    
    
    pubchem_document_term = scipy.sparse.load_npz('data/pubchem_document_term.npz'); print(pubchem_document_term.shape)
    CID_CIDidx = pd.read_pickle('CID_CIDidx.pkl'); print(CID_CIDidx.shape)
    CID_CIDidx_MeSH = pd.read_pickle('CID_CIDidx_MeSH.pkl'); print(CID_CIDidx_MeSH.shape)
    CID_CIDidx_Title = pd.read_pickle('CID_CIDidx_Title.pkl'); print(CID_CIDidx_Title.shape)
    CID_CIDidx['pubchem_idx'] = CID_CIDidx.index
    CID_CIDidx['pubchem_df'] = pubchem_document_term.sum(axis=0).A1
    CID_CIDidx['pubchem_log10df'] = np.log10(CID_CIDidx['pubchem_df'])
    CID_CIDidx['MeSH'] = CID_CIDidx_MeSH['MeSH'].str.lower()
    CID_CIDidx['Title'] = CID_CIDidx_Title['Title'].str.lower()
    #CID_CIDidx['term'] = CID_CIDidx['Title']
    CID_CIDidx = CID_CIDidx.drop(['CIDidx'], axis=1)
    CID_CIDidx.to_pickle('data/pubchem_term_info.pkl')    
    
    pubchem_document_term = scipy.sparse.load_npz('data/pubchem_document_term.npz'); print(pubchem_document_term.shape)
    PMID_PMIDidx = pd.read_pickle('PMID_PMIDidx.pkl'); print(PMID_PMIDidx.shape)
    PMID_PMIDidx = PMID_PMIDidx.drop(['PMIDidx'], axis=1)
    PMID_PMIDidx['pubchem_idx'] = PMID_PMIDidx.index
    PMID_PMIDidx['pubchem_wc'] = pubchem_document_term.sum(axis=1).A1
    PMID_PMIDidx['pubchem_log10wc'] = np.log10(PMID_PMIDidx['pubchem_wc'])
    PMID_PMIDidx.to_pickle('data/pubchem_document_info.pkl')    





def init_clinical_term_dictionary():
    patosnomed = pd.read_csv('data/patoSnoMed.txt', delimiter='\t')
    patosnomed = patosnomed[['SKSkode','Kodetekst']]
    patosnomed['SKSkode'] = patosnomed['SKSkode'].str.upper()
    patosnomed['SKSkode_rstrip'] = patosnomed['SKSkode'].str.rstrip('0')
    patosnomed['SKSkode_rstrip_len'] = patosnomed['SKSkode_rstrip'].str.len()
    patosnomed['axis'] = patosnomed['SKSkode'].str[0]
    patosnomed.Kodetekst = patosnomed.Kodetekst.str.lower()
    patosnomed['wc'] = [len(i) for i in patosnomed.Kodetekst.str.split()]
    patosnomed.to_csv('data/patosnomed.csv.gz', compression='gzip', index=None,)

    count_vectorizer = CountVectorizer(lowercase=True)
    count_vectorizer.fit_transform(patosnomed)
    pd.DataFrame(count_vectorizer.get_feature_names(), columns=['word']).to_csv('data/patosnomed_words.csv.gz', compression='gzip', index=None)
    patosnomed_words = pd.read_csv('data/patosnomed_words.csv.gz')['word'].tolist()
    patosnomed_words[0:5]    

    sct = pd.read_csv('data/sct2_Description_Full-da_DK1000005_20180930.txt', delimiter='\t')
    sct = sct[['term']]
    sct.term = sct.term.str.lower()
    sct['len'] = [len(i) for i in sct.term.str.split()]
    sct.to_csv('data/sct.csv.gz', compression='gzip', index=None,)
    sct = pd.read_csv('data/sct.csv.gz')['term'].tolist()
    sct[0:5]    
    
    count_vectorizer = CountVectorizer(lowercase=True)
    count_vectorizer.fit_transform(sct)
    pd.DataFrame(count_vectorizer.get_feature_names(), columns=['word']).to_csv('data/sct_words.csv.gz', compression='gzip', index=None)
    sct_words = pd.read_csv('data/sct_words.csv.gz')['word'].tolist()
    sct_words[0:5]    





def init_clinical_term_translation():
    sct = pd.read_csv('data/SnomedCT_InternationalRF2_PRODUCTION_20210731T120000Z/Full/Terminology/sct2_Description_Full-en_INT_20210731.txt', sep='\t'); print(sct.shape)
    sct.to_pickle('data/sct_en_term_info.pkl'); print(sct.shape)
    sct = pd.read_pickle('data/sct_en_term_info.pkl'); print(sct.shape)
    sct = sct[sct['active']==1]
    sct = sct[sct['typeId']==900000000000013009]
    sct = sct[sct['term'].str.isnumeric()==False]
    sct['term'] = sct['term'].str.strip()
    sct['term'] = sct['term'].str.replace('     ', ' ')
    sct['term'] = sct['term'].str.replace('    ', ' ')
    sct['term'] = sct['term'].str.replace('   ', ' ')
    sct['term'] = sct['term'].str.replace('  ', ' ')
    sct['term'] = sct['term'].str.lower()
    sct['ngram'] = sct['term'].str.split().apply(len)
    sct = sct.drop_duplicates(subset=['conceptId','term'])
    sct = sct[sct['ngram']<=4]
    sct.to_pickle('data/sct_en_term_info.pkl'); print(sct.shape)
    
    sct = pd.read_csv('data/SnomedCT_ManagedServiceDK_PRODUCTION_DK1000005_20210930T120000Z/Full/Terminology/sct2_Description_Full-da_DK1000005_20210930.txt', sep='\t'); print(sct.shape)
    sct.to_pickle('data/sct_da_term_info.pkl')
    sct = pd.read_pickle('data/sct_da_term_info.pkl'); print(sct.shape)
    sct = sct[sct['active']==1]
    sct = sct[sct['typeId']==900000000000013009]
    sct = sct[sct['term'].str.isnumeric()==False]
    sct['term'] = sct['term'].str.strip()
    sct['term'] = sct['term'].str.replace('    ', ' ')
    sct['term'] = sct['term'].str.replace('   ', ' ')
    sct['term'] = sct['term'].str.replace('  ', ' ')
    sct['term'] = sct['term'].str.lower()
    sct['ngram'] = sct['term'].str.split().apply(len)
    sct = sct.drop_duplicates(subset=['conceptId','term'])
    sct = sct[sct['ngram']<=4]
    sct.to_pickle('data/sct_da_term_info.pkl')    





def proximity_score():
    snomedct_term_info = pd.read_pickle('data/snomedct_term_info.pkl'); print(snomedct_term_info.shape)
    snomedct_term_info = snomedct_term_info.drop_duplicates(subset='term')
    snomedct_term_info = snomedct_term_info.reset_index(drop=True)
    snomedct_term_info = snomedct_term_info[['code', 'term']]

    snomedct_en_term_info = pd.read_pickle('data/snomedct_en_term_info.pkl'); print(snomedct_en_term_info.shape)
    snomedct_en_term_info = snomedct_en_term_info.drop_duplicates(subset='code')
    snomedct_en_term_info = snomedct_en_term_info.reset_index(drop=True)
    snomedct_en_term_info = snomedct_en_term_info[['code', 'term']]
    
    lung_records_terms = pd.read_csv('data/lung_records_terms.tsv', sep='\t', index_col=0)[['term']]; print(lung_records_terms.shape)
    
    lung_records_terms = lung_records_terms.merge(snomedct_term_info, on='term', how='left')
    lung_records_terms = lung_records_terms.merge(snomedct_en_term_info, on='code', how='left')
    
    lung_records_terms.columns = ['term_da', 'code', 'term_en']
    lung_records_terms['term_da'][lung_records_terms['code'].isna()==True]=np.nan
    lung_records_topic_modelling = pd.read_csv('data/lung_records_topic_modelling.tsv', sep='\t', index_col=0); print(lung_records_topic_modelling.shape)
    lung_records_terms['topic'] = lung_records_topic_modelling.to_numpy().argmax(axis=0)
    lung_records_terms = lung_records_terms[lung_records_terms['topic']==20]
    lung_records_terms = lung_records_terms[lung_records_terms['term_en'].isna()==False]
    lung_aging_topic_20_terms = lung_records_terms['term_en'].tolist()
    lung_aging_topic_20_terms.append('lung') #patoSnoMed:T28000:Lunge
    lung_aging_topic_20_terms

    pubmed_document_term = scipy.sparse.load_npz('data/pubmed_sct_en_document_term.npz'); print(pubmed_document_term.shape)
    pubmed_document_info = pd.read_pickle('data/pubmed_document_info.pkl'); print(pubmed_document_info.shape)
    pubmed_term_info = pd.read_pickle('data/pubmed_sct_en_term_info.pkl'); print(pubmed_term_info.shape)
    pubchem_term_info = pd.read_pickle('data/pubchem_term_info.pkl'); print(pubchem_term_info.shape)
    
    pubmed_term_info_res1 = pubmed_term_info.copy()
    pubmed_term_info_res1 = pubmed_term_info_res1[pubmed_term_info['term'].isin(pubchem_term_info['MeSH'].tolist())]
    pubmed_term_info_res1 = pubmed_term_info_res1[pubmed_term_info_res1['pubmed_df']>9]
    pubmed_term_info_res1 = pubmed_term_info_res1.reset_index(drop=True)
    
    pubmed_term_info_res2 = pubmed_term_info.copy()
    pubmed_term_info_res2 = pubmed_term_info_res2[pubmed_term_info_res2['term'].isin(lung_aging_topic_20_terms)]
    pubmed_term_info_res2 = pubmed_term_info_res2.reset_index(drop=True)
    
    filter_year = pubmed_document_info[pubmed_document_info['year']>2000]['pubmed_idx'].to_numpy()
    pubmed_document_term_filter = pubmed_document_term[filter_year,:]   
    
    pm1 = pubmed_document_term_filter[:,pubmed_term_info_res1['pubmed_idx'].to_numpy()]
    pm2 = pubmed_document_term_filter[:,pubmed_term_info_res2['pubmed_idx'].to_numpy()]
    pm1tfidf = TfidfTransformer(norm=None).fit_transform(pm1)
    pm2tfidf = TfidfTransformer(norm=None).fit_transform(pm2)
    dist_tfidf = cosine_distances(pm1tfidf.transpose(), pm2tfidf.transpose())
    dist_tfidf[dist_tfidf==1]=0

    pubmed_term_info_res = pubmed_term_info_res1.copy()
    pubmed_term_info_res['dist_tfidf'] = np.array([dist_tfidf[i,:][np.nonzero(dist_tfidf[i,:])[0]].mean() for i in range(dist_tfidf.shape[0])])
    pubmed_term_info_res['score'] = pubmed_term_info_res['dist_tfidf']
    pubmed_term_info_res = pubmed_term_info_res.sort_values(by='score')
    pubmed_term_info_res = pubmed_term_info_res.reset_index(drop=True)
    pubmed_term_info_res.to_csv('scored_terms.csv') 


# # Run code




metrics = ['euclidean', 'cosine', 'manhattan', 'l1', 'l2', 'cityblock']
matrices = ['onehot', 'tfidf', 'compressed_tfidf', 'lda_components']
methods = ['tsne', 'umap']


# # Distance matrices




for matrix in matrices:
    for metric in metrics:
        run_distance_matrix(tissue=tissue, matrix=matrix, metric=metric)


# # Embeddings




for matrix in matrices:
    for metric in metrics:
        run_umap_embedding(tissue=tissue, matrix=matrix, metric=metric)
        run_tsne_embedding(tissue=tissue, matrix=matrix, metric=metric)





for matrix in matrices:
    for metric in metrics:
        for method in methods:
            plot_embedding_mean_age(tissue=tissue, matrix=matrix, metric=metric, method=method)
            plot_embedding_cox_coef(tissue=tissue, matrix=matrix, metric=metric, method=method)
            plot_embedding_lda_topic(tissue=tissue, matrix=matrix, metric=metric, method=method)                           


# # PCA




run_pca(tissue=tissue, matrix='compressed_tfidf')
plot_pca(tissue=tissue, age_start=0, age_end=100)


# # Regression




run_linear_regression(tissue=tissue, matrix='onehot')
run_mpl_regression(tissue=tissue, matrix='onehot', max_iter=10)
run_mpl_regression(tissue=tissue, matrix='onehot', max_iter=100)
run_mpl_regression(tissue=tissue, matrix='onehot', max_iter=1000)





plot_regression(tissue=tissue, model_type='LinearRegression')





plot_regression(tissue=tissue, model_type='LinearRegression')
plot_regression_per_topic(tissue=tissue, model_type='LinearRegression')
for i in [10, 100, 1000]:
    plot_regression(tissue=tissue, model_type='MLPRegressor'+str(i))
    plot_regression_per_topic(tissue=tissue, model_type='MLPRegressor'+str(i))


# # PFI




run_permutation_importance(tissue=tissue, model_type='LinearRegression', max_samples=10000)
run_permutation_importance_per_topic(tissue=tissue, model_type='LinearRegression', n_repeats=5, n_clust=60, n_samples=10000)
for i in [10, 100, 1000]:
    run_permutation_importance(tissue=tissue, model_type='MLPRegressor'+str(i), max_samples=10000)
    run_permutation_importance_per_topic(tissue=tissue, model_type='MLPRegressor'+str(i), n_repeats=5, n_clust=60, n_samples=10000)





plot_permutation_importance(tissue=tissue, model_type='LinearRegression')
plot_permutation_importance_per_topic(tissue=tissue, model_type='LinearRegression',topic=20)
for i in [10, 100, 1000]:
    plot_permutation_importance(tissue=tissue, model_type='MLPRegressor'+str(i))
    plot_permutation_importance_per_topic(tissue=tissue, model_type='MLPRegressor'+str(i),topic=20)


# # Topic slopes




run_topic_slopes(tissue=tissue, model_type='LinearRegression')
plot_topic_slopes(tissue=tissue, model_type='LinearRegression', topic=20)
for i in [10, 100, 1000]:
    run_topic_slopes(tissue=tissue, model_type='MLPRegressor'+str(i))
    plot_topic_slopes(tissue=tissue, model_type='MLPRegressor'+str(i), topic=20)


# # Topic Cox




plot_cox(topic=20)


# # Pubmed cross referencing




init_clinical_term_translation()
init_pubmed()
init_pubchem()
proximity_score()

