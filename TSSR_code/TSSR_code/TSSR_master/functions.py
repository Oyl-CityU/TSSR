#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import os
import numpy as np
import csv
import scipy.linalg as LA
from collections import defaultdict


def load_data_from_file(dataset, folder):
    with open(os.path.join(folder, dataset+"_admat_dgc.csv"), "r") as inf: # the lncrna-disease similarity file
        read = csv.reader(inf)
        next(read)
        int_array = []
        for line in read:
            int_array.append(line[1:])

    with open(os.path.join(folder, dataset+"_simmat_dc.csv"), "r") as inf:  # the disease similarity file
        disease_sim = []
        read = csv.reader(inf)
        next(read)
        for line in read:
            disease_sim.append(line[1:])

    with open(os.path.join(folder, dataset+"_simmat_dg.csv"), "r") as inf:  # the lncrna similarity file
        lncrna_sim = []
        read = csv.reader(inf)
        next(read)
        for line in read:
            lncrna_sim.append(line[1:])


    intMat = np.array(int_array, dtype=np.double).T    # lncrna-disease interaction matrix
    diseaseMat0 = np.array(disease_sim, dtype=np.double)      #  disease semantic similarity matrix
    lncrnaMat0 = np.array(lncrna_sim, dtype=np.double)  # lncrna functional  similarity matrix

    return intMat, diseaseMat0,lncrnaMat0

def get_diseases_lncrnas_names(dataset,folder):
    with open(os.path.join(folder,dataset + "_admat_dgc.csv"),"r") as inf:
        read = csv.reader(inf)
        diseases=next(read)
        diseases=diseases[1:]
        lncenas=[line[0]for line in read]
    return diseases, lncenas

def write_metric_vector_to_file(auc_vec, file_name):
    np.savetxt(file_name, auc_vec, fmt='%.6f')


def load_metric_vector(file_name):
    return np.loadtxt(file_name, dtype=np.float64)

def gKernel(intMat):
    nl, nd = intMat.shape
    pkl = np.zeros((nl, nl))
    pkd = np.zeros((nd, nd))

    # lncRNA
    sl = np.zeros(nl)
    for i in np.arange(nl):
        sl[i] = LA.norm(intMat[i, :], 2) ** 2
    gamal = np.sum(sl) / nl
    #gamal=nl/np.sum(sl)
    for i in np.arange(nl):
        for j in np.arange(nl):
            pkl[i, j] = np.exp(-gamal * LA.norm(intMat[i, :] - intMat[j, :]) ** 2)
    # disease (这个没有用到)
    sd = np.zeros(nd)
    for i in np.arange(nd):
        sd[i] = LA.norm(intMat[:, i], 2) ** 2
    #gamad=nd/np.sum(sd)
    gamad=np.sum(sd)/nd
    for i in np.arange(nd):
        for j in np.arange(nd):
            pkd[i, j] = np.exp(-gamad * LA.norm(intMat[:, i] - intMat[:, j]) ** 2)
    return pkl, pkd
