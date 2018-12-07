def patch_desitarget_bits(desitarget_module,new_yaml):
    """
    Loads a different version of targetmask_yaml into the 
    live desitarget module passed as the first argument.
    """
    import yaml
    
    # This is a funny kind of import, not strictly necessary
    # i.e. from desitarget.targetmask import _load_mask_priorities, BitMask
    _load_mask_priorities = desitarget_module.targetmask._load_mask_priorities
    BitMask = desitarget_module.targetmask.BitMask
    
    prename = ''
    
    with open(new_yaml) as fx:
        bitdefs = yaml.load(fx)
        try:
            bitdefs = _load_mask_priorities(bitdefs,handle="priorities",prename=prename)
        except TypeError:
            pass
        try:
            bitdefs = _load_mask_priorities(bitdefs,handle="numobs",prename=prename)
        except TypeError:
            pass
        
    _bitdefs = bitdefs
    try:
        desi_mask     = BitMask('desi_mask', _bitdefs)
        mws_mask      = BitMask('mws_mask',  _bitdefs)
        bgs_mask      = BitMask('bgs_mask',  _bitdefs)
        obsconditions = BitMask('obsconditions', _bitdefs)
        obsmask       = BitMask('obsmask',       _bitdefs)
        targetid_mask = BitMask('targetid_mask', _bitdefs)
    except TypeError:
        desi_mask     = object()
        mws_mask      = object()
        bgs_mask      = object()
        obsconditions = object()
        obsmask       = object()
        targetid_mask = object()
        
    desitarget_module.targetmask.desi_mask     = desi_mask
    desitarget_module.targetmask.bgs_mask      = bgs_mask
    desitarget_module.targetmask.mws_mask      = mws_mask
    desitarget_module.targetmask.obsmask       = obsmask
    desitarget_module.targetmask.obsconditions = obsconditions
    
    return desi_mask, bgs_mask, mws_mask, obsmask, obsconditions