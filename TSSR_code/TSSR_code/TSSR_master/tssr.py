#!/usr/bin/python
#coding:utf-8
from numpy import *
import numpy as np
from scipy.stats import *

import os
import scipy.linalg as LA

class TSSR:

    def __init__(self, lambda_d=2**(-5), lambda_t=2**(-5), beta= 2**(2),max_iter=150):
        self.lambda_d = lambda_d  # default 100
        self.lambda_t = lambda_t
        self.beta=beta
        self.max_iter = max_iter

    def fix_model(self, W, intMat, drugMat, targetMat):
        size = intMat.shape
        Y=np.multiply(W,intMat)
        Sd=drugMat
        St=targetMat
        lowest_score=Inf
        n = size[0]
        m = size[1]
        A_old = drugMat
        B_old = targetMat
        yy = (A_old.sum(axis=1))
        yy1 = yy.reshape(n, 1)
        zz = (B_old.sum(axis=0))
        zz1 = zz.reshape(1, m)
        A_old = A_old / np.tile(yy1, (1, n))
        B_old = B_old / np.tile(zz1, (m, 1) )
        for i in xrange(self.max_iter):
            B_old = np.mat(B_old)
            A_old = np.mat(A_old)
            Y,W=np.mat(Y),np.mat(W)
            Sd,St=np.mat(Sd),np.mat(St)
            Dan=2*np.multiply(W,Y)*(B_old.T*Y.T)+2*self.lambda_d*Sd
            Dap=2*(np.multiply(W,A_old*Y*B_old))*B_old.T * Y.T+2*self.lambda_d * A_old + self.beta*np.ones([n,n])
            Dapeps=np.finfo(float).eps * ones([n,n])
            cc=np.array((A_old / (Dap +Dapeps)).sum(axis=1).T)
            aa = np.diag(cc[0])
            ba=(np.multiply(A_old,Dan/(Dap+Dapeps ))).sum(axis=1)
            A=np.multiply(A_old,aa*Dan + np.ones([n,n]))/(aa*Dap + np.tile(ba,(1,n))+Dapeps)
            Dbn=2*(Y.T * A.T)* (np.multiply(W,Y))+ 2 * self.lambda_t * St
            Dbp=2 * Y.T * A.T * (np.multiply(W,A*Y*B_old))+ 2* self.lambda_t * B_old + self.beta*np.ones([m,m])
            Dbpeps=np.finfo(float).eps * ones([m,m])
            dd=np.array((B_old/(Dbp+ Dbpeps)).sum(axis=0))
            ab=np.diag(dd[0])
            bb=(np.multiply(B_old,Dbn/(Dbp+Dbpeps))).sum(axis=0)
            B= np.multiply(B_old,(Dbn*ab +np.ones([m,m]))/(Dbp* ab + np.tile(bb,(m,1))+ Dbpeps))

            if (abs(A-A_old)).sum(axis=0).sum(axis=1)+ (abs(B-B_old)).sum(axis=0).sum(axis=1)< 1e-6:
                break
            A_old=A
            B_old=B

        score=np.linalg.norm(np.multiply(W,Y-A_old * Y * B_old),'fro')**2 + self.lambda_d * np.linalg.norm(Sd-A_old,'fro')**2 +  self.lambda_t * np.linalg.norm(St-B_old,'fro')**2
        + self.beta * (A_old.sum(axis=0).sum(axis=1))+ self.beta * (B_old.sum(axis=0).sum(axis=1))

        if score < lowest_score:
            A_final=(A+A.T)/2
            B_final=(B+B.T)/2
            lowest_score=score

        Y_pre= A_final * Y *  B_final
        self.predictR= np.array(Y_pre)


    def predict_scores(self, test_data, N):
        inx = np.array(test_data)
        return self.predictR[inx[:, 0], inx[:, 1]]


    def __str__(self):
        return "Model: TSSR,  lambda_d:%s, lambda_t:%s, beta:%s, max_iter:%s" % ( self.lambda_d, self.lambda_t,self.beta, self.max_iter)