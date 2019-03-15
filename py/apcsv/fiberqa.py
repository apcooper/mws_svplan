import sys
import os
import numpy as np
import glob
from   astropy.table import Table
from   astropy.io    import fits

import matplotlib.pyplot as pl

import apcsv.util as util

from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

############################################################
def glob_fibers(path):
    """
    """
    return glob.glob(os.path.join(path,'tile-*.fits'))

############################################################
def load_fiber_data(path,ext='FIBERASSIGN'):
    """
    """
    data = dict()
    fiber_files = glob_fibers(path)
    print('Have {} tiles'.format(len(fiber_files)))
    for fiber_file in fiber_files:
        filename = os.path.basename(fiber_file)
        itile    = int(os.path.splitext(filename)[0].split('-')[1])
        d        = fits.getdata(fiber_file,ext)
        data[itile] = d
    return data

############################################################
def fibers_assigned(tiledata):
    """
    """
    # Also 'FIBERSTATUS', 'FIBERMASK'
    return tiledata['DESI_TARGET'] > 0

############################################################
def format_cats(cats,ncols=60):
    """
    """
    fmt_string = '{{:{ncols:d}}}'.format(ncols=ncols)
    cats    = ' & '.join(cats)
    cat_str = list()
    while len(cats) > ncols:
        rindex = cats[0:ncols].rindex(' & ')+1
        cat_str.append(fmt_string.format(cats[0:ncols][0:rindex]))
        cats = cats[rindex:]
    cat_str.append(fmt_string.format(cats))
    
    return cat_str
    
############################################################
def report_fibers_mws(tiledata,ncols=60):
    """
    """
    fmt_string   = '{{:{ncols:d}}} {{:<10d}}'.format(ncols=ncols)
    non_mws_line = None
    
    lines = list()
    ubits, ubits_count = np.unique(tiledata['MWS_TARGET'],return_counts=True)
    for ubit, ubit_count in zip(ubits, ubits_count):
        cats = list()
        for _ in mws_mask.names(ubit):
            _ = _.replace('MWS_','')
            cats.append(_)

        cat_str = format_cats(cats)
        if len(cat_str) > 1:
            for _s in cat_str[:-1]:
                lines.append(_s)
                
        cat_str = cat_str[-1]
        blank   = len(cat_str.strip()) == 0
        
        if blank: cat_str = '(non-MWS targets)'
            
        line = fmt_string.format(cat_str,ubit_count)
        
        if blank:
            non_mws_line = line
        else:
            lines.append(line)
    
    for line in lines:
        print(line)
    
    if non_mws_line is not None:
        print()
        print(non_mws_line)

    return

############################################################
def plot_pie_target_classes(tiledata):
    """
    """
    ubits, ubits_count = np.unique(tiledata['MWS_TARGET'],return_counts=True)
    
    labels = [ubits_labels[u] if u in util.mws_bits_to_name else u for u in ubits]
    colors = [ubits_colors[u] if u in util.mws_bits_to_colors else 'None' for u in ubits]

    kwargs = dict()
    kwargs['labels'] = labels
    kwargs['colors'] = colors

    pl.pie(ubits_count,**kwargs)
    pl.axis('equal')
    
    return

############################################################
def survey_overlaps(tiledata):
    """
    Returns True for fibers that have MWS_ANY bit set in 
    DESI_TARGET and other DESI_TARGET bits set.
    """
    is_assigned = fibers_assigned(tiledata)
    
    # Define as having MWS_ANY and any bits other than MWS_ANY
    is_mws   = (tiledata['DESI_TARGET'] &  desi_mask.mask('MWS_ANY'))!= 0
    is_other = (tiledata['DESI_TARGET'] & ~desi_mask.mask('MWS_ANY'))!= 0
    return (is_assigned & is_mws & is_other)

############################################################
def report_survey_overlaps(tiledata):
    """
    """
    is_overlap = survey_overlaps(tiledata)
    print('{:d} fibers with both MWS and non-MWS target bits:'.format(is_overlap.sum()))
    
    udesi,udesi_count = np.unique(tiledata['DESI_TARGET'][is_overlap],
                                       return_counts=True)
    
    overlap_lines = list()
    for (_desi,_c) in zip(udesi,udesi_count):
        names = desi_mask.names(_desi)
        names.remove('MWS_ANY')
        desi_str = ' & '.join(names)
        
        _ = tiledata['DESI_TARGET'][is_overlap] == _desi
        mws_bits = np.unique(tiledata['MWS_TARGET'][is_overlap][_])
        
        mws_str = list()
        for _ in mws_bits:
            if _ in util.mws_bits_to_name:
                mws_str.append(util.mws_bits_to_name[_])
        mws_str = ' | '.join(mws_str)
        
        line = ' {:30s} {:4d} {:s}'.format(desi_str,_c,mws_str)
        overlap_lines.append(line)
        
    for line in (sorted(overlap_lines)):
        print(line)
        
    return