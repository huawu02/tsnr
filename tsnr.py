#!/usr/bin/env python

import os
import nibabel as nb
import numpy as np
#from nipype.interfaces import fsl

if __name__ == '__main__':

    import argparse

    arg_parser = argparse.ArgumentParser()
    arg_parser.description  = ('Calculate temporal SNR of a 4D nifti that contains a time series.\n\n')
    arg_parser.add_argument('infile', help='path to nifti file of the time series')
    arg_parser.add_argument('-o', '--outbase', default='', help='basename of the output files')
    arg_parser.add_argument('-d', '--discard_vol', type=int, default=3, help='number of volumes to discard from the beginning')
    arg_parser.add_argument('-f', '--bet_frac', type=float, default=0.5, help='fraction for FSL''s BET mask generation')
    args = arg_parser.parse_args()

    basename = (os.path.basename(args.infile)).split('.')[0]
    outbase = args.outbase
    if args.outbase == '':
        outbase = basename+'_tsnr'
    masked_name = basename+'_masked.nii.gz'
    mask_name = basename+'_masked_mask.nii.gz'

    #fsl.BET(in_file=args.infile, frac=args.bet_frac, mask=True, no_output=False, out_file=masked_name).run()
    os.system("fsl5.0-bet %s %s -f %f -m" %(args.infile, masked_name, args.bet_frac)) 
    mask = nb.load(mask_name).get_data() > 0.5

    ni = nb.load(args.infile)
    d = ni.get_data()
    d = d[...,args.discard_vol:]  # discard the first volumes
    tsnr = np.multiply(np.mean(d,axis=3) / np.std(d,axis=3), mask)
    tsnr[np.isnan(tsnr)] = 0
    ni_snr = nb.Nifti1Image(tsnr, ni.get_affine())
    nb.save(ni_snr, outbase+'.nii.gz')

    os.remove(masked_name)
    os.remove(mask_name)
