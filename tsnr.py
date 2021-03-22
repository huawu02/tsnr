#!/usr/bin/env python

# Compute temperal SNR of a time series, with the option to remove low frequency temporal fluctuations
# and compute SFNR (signal fluctuation to noise ratio). Also perform Weisskoff analysis and calculate 
# the STD within ROIs.
# Require AFNI functions 3dTcat, 3dAutomask, 3dDetrend.

import os
import json
import nibabel as nb
import numpy as np
from scipy import ndimage
# disable X-display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# from nipype.interfaces import fsl

if __name__ == '__main__':

    import argparse

    arg_parser = argparse.ArgumentParser()
    arg_parser.description  = ('Calculate temporal SNR of a 3D+time nifti after removing components from voxel time series.\n\n')
    arg_parser.add_argument('infile', help='path to nifti file of the time series')
    arg_parser.add_argument('-o', '--outbase', default='', help='basename of the output files')
    arg_parser.add_argument('-d', '--discard_vol', type=int, default=3, help='number of volumes to discard from the beginning, default=3')
    arg_parser.add_argument('-f', '--mask_frac', type=float, default=0.4, help='clip level fraction for mask generation, default=0.4')
    arg_parser.add_argument('-p', '--detrend_polort', type=int, default=2, help='polynomials order for 3dDetrend, default=2')
    arg_parser.add_argument('-r', '--roi_size', type=int, default=21, help='length of square ROI in Weisskoft analysis, default=21')
    arg_parser.add_argument('--save_all_outputs', action='store_true', help='flag for saving all intermediate results')
    args = arg_parser.parse_args()

    basename = (os.path.basename(args.infile)).split('.')[0]
    outbase = args.outbase
    if outbase == '':
        outbase = basename
    reorient_name = outbase+'_reorient.nii.gz'
    tseries_name = outbase+'_tseries.nii.gz'
    detrend_name = outbase+'_detrend.nii.gz'
    mask_name = outbase+'_mask.nii.gz'
    tsnr_name = outbase+'_tsnr.nii.gz'
    tmean_name = outbase+'_tmean.nii.gz'
    tstd_name = outbase+'_tstd.nii.gz'
    nstd_name = outbase+'_nstd.nii.gz'
    sfnr_name = outbase+'_sfnr.nii.gz'

    # reorient the infile to standard RAS coordinates
    # fsl.Reorient2Std(in_file=args.infile, out_file=reorient_name).run()
    cmd = "fsl5.0-fslreorient2std " + "%s %s" % (args.infile, reorient_name)
    status = os.system(cmd)
    
    
    # discard volumes, mask, detrend
    # use the first line when runnig in local env, use second line when running in docker
    # os.system("3dTcat -prefix %s %s[%d..$]; 3dAutomask -prefix %s -clfrac %f %s; 3dDetrend -prefix %s -polort %d %s" 
    os.system(". /etc/afni/afni.sh; 3dTcat -prefix %s %s[%d..$]; 3dAutomask -prefix %s -clfrac %f %s; 3dDetrend -prefix %s -polort %d %s" 
              %(tseries_name, reorient_name, args.discard_vol, 
                mask_name, args.mask_frac, tseries_name, 
                detrend_name, args.detrend_polort, tseries_name))

    ni = nb.load(tseries_name)
    pixdim = ni.get_header()['pixdim'][1:4]  # pixel size in mm
    tr = ni.get_header()['pixdim'][4]  # TR in seconds
    affine = ni.affine
    tseries = np.asarray(ni.dataobj)
    noise = nb.load(detrend_name).get_data()
    mask  = nb.load(mask_name).get_data()
    tmean = np.multiply(np.mean(tseries, axis=3), mask)
    tstd  = np.multiply(np.std(tseries, axis=3), mask)
    nstd  = np.multiply(np.std(noise, axis=3), mask)
    tsnr  = np.nan_to_num(np.multiply(np.divide(tmean, tstd), mask))
    sfnr  = np.nan_to_num(np.multiply(np.divide(tmean, nstd), mask))
    
    # center of mass drift 
    num_tpoints = noise.shape[3]
    center_of_mass_tseries_coord = np.zeros((3, num_tpoints))
    center_of_mass_tseries_pixel = np.zeros((3, num_tpoints))
    center_of_mass_drift = np.zeros((3, num_tpoints))
    for t in range(num_tpoints):
        center_of_mass_tseries_pixel[:, t] = ndimage.measurements.center_of_mass(tseries[...,t])     # pixel index of the center of mass
        center_of_mass_tseries_coord[:, t] = np.dot(affine[0:3, 0:3], center_of_mass_tseries_pixel[:, t].T) + affine[3, 0:3]   # coordinates of the center of mass
        if t == 0:
            center_of_mass = list(int(s) for s in ndimage.measurements.center_of_mass(tseries[...,t]))
        else:
            center_of_mass_drift[:, t] = center_of_mass_tseries_coord[:, t] - center_of_mass_tseries_coord[:, 0]
    t = np.arange(num_tpoints) * tr
    t = t[:, np.newaxis]
    drift = np.zeros((num_tpoints, 1))
    drift_rate = np.zeros((3, 1))
    for i in range(3):
        drift[:,0] = center_of_mass_drift[i,:].T
        drift_rate[i], _, _, _ = np.linalg.lstsq(t, drift)   # rate of drift in mm

    # Weisskoff analysis radius of decorrelation
    roi_length = range(1, args.roi_size+1)
    roi_std_detrend = []
    for r in roi_length:
        roi = []
        for i in range(0,2):
            roi.append(range(center_of_mass[i]-r//2, center_of_mass[i]+r//2+np.mod(r,2)))
        roi_mask = np.zeros(noise.shape)
        roi_mask[np.meshgrid(roi[0],roi[1],[center_of_mass[2]],range(num_tpoints))] = 1
        #roi_std = np.std(tseries[np.where(roi_mask)].flatten())
        #roi_std_noise = np.std(np.sum(np.sum(np.multiply(noise[:,:,center_of_mass[2],:], roi_mask[:,:,center_of_mass[2],:]), axis=0), axis=0)*1.0/r/r)
        roi_mean = np.sum(np.sum(np.multiply(tseries[:,:,center_of_mass[2],:], roi_mask[:,:,center_of_mass[2],:]), axis=0), axis=0)*1.0/r/r
        coeff = np.polyfit(range(num_tpoints), roi_mean, 2)
        poly = np.poly1d(coeff)
        res = roi_mean - poly(range(num_tpoints))
        roi_std_detrend.append(np.std(res))
    rdc = roi_std_detrend[0] / roi_std_detrend[args.roi_size-1]
    roi_signal_mean = roi_mean / np.mean(roi_mean) # normalize to the temporal mean signal
    roi_signal_mean_fitted = poly(range(num_tpoints)) / np.mean(roi_mean)
    roi_residual = roi_signal_mean - roi_signal_mean_fitted
    roi_temporal_variance = np.std(roi_residual) 
    sfnr_center = np.mean(sfnr[np.where(roi_mask[:,:,:,0])])
    sfnr_edge = np.percentile(sfnr[:,:,center_of_mass[2]][np.where(mask[:,:,center_of_mass[2]])], 95)
    tsnr_center = np.mean(tsnr[np.where(roi_mask[:,:,:,0])])
    tsnr_edge = np.percentile(tsnr[:,:,center_of_mass[2]][np.where(mask[:,:,center_of_mass[2]])], 95)
    
    # save data
    data = {'roi_size': args.roi_size, 
            'roi_std': ['%.4f' % x for x in roi_std_detrend], 
            'radius_decorrelation': '%.1f' % rdc,
            'roi_signal_mean': ['%.4f' % x for x in roi_signal_mean],
            'roi_signal_mean_fitted': ['%.4f' % x for x in roi_signal_mean_fitted],
            'roi_residual': ['%.4f' % x for x in roi_residual],
            'roi_temporal_variance': '%.6f' % roi_temporal_variance,
            'center_of_mass_x': ['%.4f' % x for x in center_of_mass_tseries_coord[0, :]],
            'center_of_mass_y': ['%.4f' % x for x in center_of_mass_tseries_coord[1, :]],
            'center_of_mass_z': ['%.4f' % x for x in center_of_mass_tseries_coord[2, :]],
            'center_of_mass_drift_x': ['%.4f' % x for x in center_of_mass_drift[0, :]],
            'center_of_mass_drift_y': ['%.4f' % x for x in center_of_mass_drift[1, :]],
            'center_of_mass_drift_z': ['%.4f' % x for x in center_of_mass_drift[2, :]],
            'center_of_mass_drift_rate': ['%.6f' % x for x in drift_rate],
            'sfnr_center': '%.2f' % sfnr_center,
            'sfnr_edge': '%.2f' % sfnr_edge,
            'tsnr_center': '%.2f' % tsnr_center,
            'tsnr_edge': '%.2f' % tsnr_edge}
    with open(outbase+'_results.json','w') as fp:
        json.dump(data, fp, indent=2, sort_keys=True)

    # plot Weisskoff analysis result
    plt.plot(range(1, args.roi_size+1), roi_std_detrend, '-x', range(1, args.roi_size+1), np.divide(roi_std_detrend[0], range(1, args.roi_size+1)), 'r-x')
    plt.axhline(roi_std_detrend[args.roi_size-1], color='k', linestyle='--', linewidth=0.5)
    plt.axvline(rdc, color='k', linestyle='--', linewidth=0.5) 
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0, args.roi_size+1)
    ymin = 0.1 if roi_std_detrend[0]/args.roi_size < 1 else 1
    plt.ylim(ymin, roi_std_detrend[0]*2)
    plt.xticks([1, 10, (args.roi_size//10)*10], [1, 10, (args.roi_size//10)*10])
    plt.yticks([0.1, 1, 10] if ymin==0.1 else [1, 10], [0.1, 1, 10] if ymin==0.1 else [1, 10])
    plt.xlabel('ROI length (N)')
    plt.ylabel('Standard deviation within ROI')
    plt.text(1.5, 1.5*roi_std_detrend[0], 'radius of decorrelation = %.1f' %(rdc))
    plt.savefig(outbase+'_rdc.png', bbox_inches='tight')
    plt.close()

    # plot ROI mean signal of the time series
    fig, ax1 = plt.subplots()
    p1, = ax1.plot((np.arange(num_tpoints) + args.discard_vol + 1) * tr, roi_signal_mean, label='mean signal')
    p2, = ax1.plot((np.arange(num_tpoints) + args.discard_vol + 1) * tr, roi_signal_mean_fitted, label='fitted')
    ax1.set_xlabel('time (s)')
    ax1.set_xlim(0, (num_tpoints+1)*tr)
    ax1.set_ylabel('signal intensity (normalized)')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

    ax2 = ax1.twinx()
    p3, = ax2.plot((np.arange(num_tpoints) + args.discard_vol + 1) * tr, roi_residual, label='residual', color='r')
    ax2.set_ylabel('residual signal (normalized)')
    lines = [p1, p2, p3]
    ax2.legend(lines, [l.get_label() for l in lines], loc=1)

    plt.text(0.1, 0.95, 'temporal variation = {:.3f}%'.format(100 * roi_temporal_variance), transform=ax1.transAxes)
    plt.savefig(outbase+'_mean_signal.png', bbox_inches='tight')
    plt.close()

    # plot drift of center of mass
    t = (np.arange(num_tpoints) + args.discard_vol + 1) * tr
    plt.plot(t, center_of_mass_drift[0, :], 'b', label='x')
    plt.plot(t, center_of_mass_drift[1, :], 'g', label='y')
    plt.plot(t, center_of_mass_drift[2, :], 'r', label='z')
    plt.plot(t, (t - (args.discard_vol+1)*tr) * drift_rate[0], 'b--')
    plt.plot(t, (t - (args.discard_vol+1)*tr) * drift_rate[1], 'g--')
    plt.plot(t, (t - (args.discard_vol+1)*tr) * drift_rate[2], 'r--')

    plt.legend(loc='upper right')
    plt.ylim(-np.ceil(np.max(np.abs(center_of_mass_drift)) * 10) / 10.0, np.ceil(np.max(np.abs(center_of_mass_drift)) * 10) / 10.0)
    plt.ylabel('drift of center of mass (mm)'); plt.xlabel('time (s)'); plt.xlim(0, (num_tpoints+1)*tr)
    plt.text(0.1, 0.999999999, 'drift rate (mm/s)\nx = %.1e \ny = %.1e \nz = %.1e' %(drift_rate[0], drift_rate[1], drift_rate[2]), transform=ax1.transAxes)
    plt.savefig(outbase+'_cm_drift.png', bbox_inches='tight')
    plt.close()

    ni_sfnr  = nb.Nifti1Image(sfnr,  ni.get_affine())
    nb.save(ni_sfnr,  sfnr_name)
    if args.save_all_outputs:
        ni_tmean = nb.Nifti1Image(tmean, ni.affine)
        ni_nstd  = nb.Nifti1Image(nstd,  ni.affine)
        ni_tstd  = nb.Nifti1Image(tstd,  ni.affine)
        ni_tsnr  = nb.Nifti1Image(tsnr,  ni.affine)
        nb.save(ni_tmean, tmean_name)
        nb.save(ni_nstd,  nstd_name)
        nb.save(ni_tstd,  tstd_name)
        nb.save(ni_tsnr,  tsnr_name)
    else:
        os.remove(reorient_name)
        os.remove(tseries_name)
        os.remove(detrend_name)
        os.remove(mask_name)

