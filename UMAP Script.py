# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 20:30:47 2018

@author: wudd1
"""
def run_umap(file):
    import umap
    import pandas as pd
    data = file
    from sklearn.datasets import load_digits
    digits = load_digits()
    embedding = umap.UMAP(n_neighbors=5, min_dist=0.3, metric='correlation').fit_transform(digits.data)

    dataset = data.values

    dataset_t = dataset
#dataset_t = dataset.transpose()

    embedding = umap.UMAP(n_neighbors=10).fit_transform(dataset_t)
#embedding = umap.UMAP(n_neighbors=10, min_dist=0.3,metric='correlation').fit_transform(dataset_t)

    rownames = list(data.index)
    colnames = ["umap dim 1","umap dim 2"]

    res = pd.DataFrame(embedding, index=rownames, columns=colnames)
    res = res.join(data)
    return res.to_csv("umap_data.txt", sep = '\t', index=True, header=True)
