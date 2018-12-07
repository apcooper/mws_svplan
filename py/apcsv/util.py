from desitarget.targetmask import mws_mask
############################################################
def expand_mws(*abbrevs):
    """
    """
    expand = {
        'M'  : 'BROAD',
        'MR' : 'MAIN_RED',
        'MB' : 'MAIN_BLUE',
        'BF' : 'MAIN_BLUE_FAINT',
        'RF' : 'MAIN_RED_FAINT',
        'MN' : 'BROAD_NORTH',
        'MS' : 'BROAD_SOUTH',
        'MRN': 'MAIN_RED_NORTH',
        'MRS': 'MAIN_RED_SOUTH',
        'MBN': 'MAIN_BLUE_NORTH',
        'MBS': 'MAIN_BLUE_SOUTH',
        'BFN': 'MAIN_BLUE_FAINT_NORTH',
        'BFS': 'MAIN_BLUE_FAINT_SOUTH',
        'RFN': 'MAIN_RED_FAINT_NORTH',
        'RFS': 'MAIN_RED_FAINT_SOUTH',
    }
    expansion = list()
    for _ in abbrevs:
        if _ in expand: 
            expansion.append('MWS_'+expand[_])
        else:
            expansion.append('MWS_'+_)
    return '|'.join(expansion)

############################################################

mws_bits_to_name = {
    0: 'Non-MWS',
    mws_mask.mask(expand_mws('M','MN','MS')): 'BROAD (RED REJECT)',
    mws_mask.mask(expand_mws('MR','MRN','MRS')): 'RED',
    mws_mask.mask(expand_mws('MB','MBN','MBS')): 'BLUE',
    mws_mask.mask(expand_mws('WD')): 'WD',
    mws_mask.mask(expand_mws('NEARBY')): 'NEARBY',
    mws_mask.mask(expand_mws('WD','NEARBY')): 'WD & NEARBY',
    mws_mask.mask(expand_mws('M','MN','MS','NEARBY')): 'BROAD & NEARBY',
    mws_mask.mask(expand_mws('WD','MB','MBN','MBS')): 'BLUE & WD',
    mws_mask.mask(expand_mws('BF','BFN','BFS')): 'BLUE FAINT',
    mws_mask.mask(expand_mws('RF','RFN','RFS')): 'RED FAINT',
    mws_mask.mask(expand_mws('WD','BF','BFN','BFS')): 'BLUE FAINT & WD',
    mws_mask.mask(expand_mws('WD','RF','RFN','RFS')): 'RED FAINT & WD',
}

mws_name_to_bits = dict([v,k] for k,v in mws_bits_to_name.items())

mws_bits_to_colors = {
        0: 'lightgrey',
        mws_name_to_bits['BROAD (RED REJECT)']: 'k',
        mws_name_to_bits['RED']: 'r',
        mws_name_to_bits['BLUE']: 'b',
        mws_mask.mask(expand_mws('WD')): 'c',
        mws_mask.mask(expand_mws('NEARBY')): 'yellow',
        mws_mask.mask(expand_mws('WD','NEARBY')): 'green',
        mws_mask.mask(expand_mws('M','MN','MS','NEARBY')): 'pink',
        mws_mask.mask(expand_mws('WD','MB','MBN','MBS')): 'grey',
}