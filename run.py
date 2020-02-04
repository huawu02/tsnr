#!/usr/bin/env python

import json
import os

if __name__ == '__main__':

    # Parse all of the input arguments from the config.json file
    config_file = '/flywheel/v0/config.json'
    if not os.path.isfile(config_file):
        raise AssertionError('No Config File FOUND!')
    else:
        with open(config_file, 'r') as f:
            config = json.load(f)

    infile = config['inputs']['nifti']['location']['path']
    basename = (os.path.basename(infile)).split('.')[0]
    discard_vol = config['config']['discarded_volume']
    mask_threshold = config['config']['mask_threshold']
    roi_size = config['config']['roi_size']
    save_all_outputs = config['config']['save_all_outputs']

    # Set output name
    outdir = '/flywheel/v0/output'
    outpath = os.path.join(outdir, basename)
    if save_all_outputs:
        os.system('python /flywheel/v0/tsnr.py %s -o %s -d %d -f %f -r %d --save_all_outputs' %(infile, outpath, discard_vol, mask_threshold, roi_size))
