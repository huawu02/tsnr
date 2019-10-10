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
    if outbase == '':
        outbase = basename
    outbase_tsnr = outbase+'_tsnr'
    outbase_tmean = outbase+'_tmean'
    outbase_noise = outbase+'_noise'
    masked_name = basename+'_masked.nii.gz'
    mask_name = basename+'_masked_mask.nii.gz'

    #fsl.BET(in_file=args.infile, frac=args.bet_frac, mask=True, no_output=False, out_file=masked_name).run()
    os.system("fsl5.0-bet %s %s -f %f -m" %(args.infile, masked_name, args.bet_frac)) 
    mask = nb.load(mask_name).get_data() > 0.5

    ni = nb.load(args.infile)
    d = ni.get_data()
    d = d[...,args.discard_vol:]  # discard the first volumes
    noise = np.multiply(np.std(d,axis=3),mask)
    tmean = np.multiply(np.mean(d,axis=3),mask)
    tsnr = np.zeros(noise.shape)
    tsnr[mask] = np.divide(tmean[mask], noise[mask])
    ni_noise = nb.Nifti1Image(noise, ni.get_affine())
    nb.save(ni_noise, outbase_noise+'.nii.gz')
    ni_tmean = nb.Nifti1Image(tmean, ni.get_affine())
    nb.save(ni_tmean, outbase_tmean+'.nii.gz')
    ni_tsnr = nb.Nifti1Image(tsnr, ni.get_affine())
    nb.save(ni_tsnr, outbase_tsnr+'.nii.gz')

    os.remove(masked_name)
    os.remove(mask_name)
