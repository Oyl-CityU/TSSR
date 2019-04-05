#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import time
import getopt
from functions import *
from tssr import TSSR

print os.getcwd()


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "m:d:f:c:s:o:n:p",
                                   ["method=", "dataset=", "data-dir=", "specify-arg=", "method-options=",
                                    "predict-num=", "output-dir=", ])
    except getopt.GetoptError:
        sys.exit()

    data_dir = os.path.join(os.path.pardir, 'data')
    output_dir = os.path.join(os.path.pardir, 'output')
    sp_arg, model_settings, predict_num = 1, [], 100
    for opt, arg in opts:
        if opt == "--method":
            method = arg
        if opt == "--dataset":
            dataset = arg
        if opt == "--data-dir":
            data_dir = arg
        if opt == "--output-dir":
            output_dir = arg
        if opt == "--specify-arg":
            sp_arg = int(arg)
        if opt == "--method-options":
            model_settings = [s.split('=') for s in str(arg).split()]
        if opt == "--predict-num":
            predict_num = int(arg)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    intMat, diseaseMat0, lncrnaMat0 = load_data_from_file(dataset, os.path.join(data_dir, 'datasets'))
    disease_names, lncrna_names = get_diseases_lncrnas_names(dataset, os.path.join(data_dir, 'datasets'))
    diseaseMat1, lncrnaMat1 = gKernel(intMat)
    diseaseMat = 0.5 * (diseaseMat0 + diseaseMat1)
    lncrnaMat = 0.5 * (lncrnaMat0 + lncrnaMat1)
    if method=="tssr":
        args={'lambda_d':2**(-5),'lambda_t':2**(-5),'beta':2**(2),'max_iter':150}
    for key, val in model_settings:
        args[key] = int(val)

    if predict_num > 0:

        if method == 'tssr':
            model = TSSR(lambda_d=args['lambda_d'], lambda_t=args['lambda_t'], beta=args['beta'],
                         max_iter=args['max_iter'])
        cmd = str(model)
        print "Dataset:" + dataset + "\n" + cmd
        W = np.ones(intMat.shape)
        model.fix_model(W, intMat, diseaseMat, lncrnaMat)
        x, y = np.where(intMat == 0)
        scores = model.predict_scores(zip(x, y), 5)
        ii = np.argsort(scores)[::-1]
        predict_pairs = [(disease_names[x[i]], lncrna_names[y[i]], scores[i]) for i in ii[:predict_num]]
        fileObject = open(os.path.join(output_dir, "_".join([method, dataset, "new_dti.txt"])), 'w')
        nn = 1
        for ii in predict_pairs:
            seq = (str(nn), ii[0], ii[1], str(ii[2]))
            fileObject.write('\t'.join(seq))
            fileObject.write('\n')
            nn += 1
        fileObject.close()


if __name__ == "__main__":
    main(sys.argv[1:])
