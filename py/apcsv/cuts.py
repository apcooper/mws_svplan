import numpy as np
import warnings 

# These are custom versions of the MWS functions in desitarget.cuts.

############################################################
def isMWS_faint_colors(gflux=None, rflux=None, zflux=None, w1flux=None, w2flux=None,
                       obs_rflux=None,
                       primary=None, south=True):
    """Set of color-based cuts used by MWS target selection classes
    (see, e.g., :func:`~desitarget.cuts.isMWS_main` for parameters).
    """
    if primary is None:
        primary = np.ones_like(rflux, dtype='?')
  
    # APC no gaia and different magitude cuts for the faint sample
    red_faint  = primary.copy()
    blue_faint = primary.copy()

    # APC faint targets are 19 <= r < 22
    red_faint  &= rflux > 10**((22.5-22.0)/2.5)
    blue_faint &= rflux > 10**((22.5-22.0)/2.5)
    
    red_faint  &= rflux <= 10**((22.5-19.0)/2.5)
    blue_faint &= rflux <= 10**((22.5-19.0)/2.5)

    # APC faint targets are robs < 23
    red_faint  &= obs_rflux > 10**((22.5-23.0)/2.5)
    blue_faint &= obs_rflux > 10**((22.5-23.0)/2.5)

    # APC color cut for BLUE-FAINT
    blue_faint &= rflux < gflux * 10**(0.7/2.5) # (g-r)<0.7
    
    # APC color cut for RED-FAINT
    red_faint &= rflux >= gflux * 10**(0.7/2.5) # (g-r)>=0.7

    return red_faint, blue_faint

############################################################
def isMWS_main_colors(gflux=None, rflux=None, zflux=None, w1flux=None, w2flux=None,
                      pmra=None, pmdec=None, parallax=None, obs_rflux=None,
                      gaiagmag=None, gaiabmag=None, gaiarmag=None, gaiaaen=None,
                      primary=None, south=True):
    """Set of color-based cuts used by MWS target selection classes
    (see, e.g., :func:`~desitarget.cuts.isMWS_main` for parameters).
    """
    if primary is None:
        primary = np.ones_like(rflux, dtype='?')
    mws = primary.copy()
   
    # ADM main targets are point-like based on GAIA_ASTROMETRIC_NOISE.
    mws &= gaiaaen < 3.0

    # ADM main targets are 16 <= r < 19
    mws &= rflux > 10**((22.5-19.0)/2.5)
    mws &= rflux <= 10**((22.5-16.0)/2.5)

    # ADM main targets are robs < 20
    mws &= obs_rflux > 10**((22.5-20.0)/2.5)

    # ADM calculate the overall proper motion magnitude
    # ADM inexplicably I'm getting a Runtimewarning here for
    # ADM a few values in the sqrt, so I'm catching it
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pm = np.sqrt(pmra**2. + pmdec**2.)

    # ADM make a copy of the main bits for a red/blue split
    red  = mws.copy()
    blue = mws.copy()

    # ADM MWS-BLUE is g-r < 0.7
    blue &= rflux < gflux * 10**(0.7/2.5)                      # (g-r)<0.7
 
    # ADM MWS-RED and MWS-BROAD have g-r >= 0.7
    red &= rflux >= gflux * 10**(0.7/2.5)                      # (g-r)>=0.7
    
    # APC BROAD astrometric cuts are the complement of those for RED
    broad = red.copy()

    # ADM MWS-RED also has parallax < 1mas and proper motion < 7.
    red &= pm < 7.
    red &= parallax < 1.

    # ADM MWS-BROAD has parallax > 1mas OR proper motion > 7.
    broad &= (parallax >= 1.) | (pm >= 7.)
    
    return broad, red, blue

############################################################
def notinMWS_main_mask(gaia=None, gfracmasked=None, gnobs=None, gflux=None,
                       rfracmasked=None, rnobs=None, rflux=None,
                       objtype=None,
                       gaiadupsource=None, primary=None):
    """Standard set of masking-based cuts used by MWS target selection classes
    (see, e.g., :func:`~desitarget.cuts.isMWS_main` for parameters).
    """
    if primary is None:
        primary = np.ones_like(gaia, dtype='?')
    mws = primary.copy()

    # ADM apply the mask/logic selection for all MWS-MAIN targets
    # ADM main targets match to a Gaia source
    mws &= gaia
    mws &= (gfracmasked < 0.5) & (gflux > 0) & (gnobs > 0)
    mws &= (rfracmasked < 0.5) & (rflux > 0) & (rnobs > 0)

    mws &= ~gaiadupsource
    
    mws &= (objtype == b'PSF')

    return mws

############################################################
def notinMWS_faint_mask(gfracmasked=None, gnobs=None, gflux=None,
                        rfracmasked=None, rnobs=None, rflux=None,
                        objtype=None,
                        primary=None):
    """
    Masking cuts for faint stars in SV
    """
    if primary is None:
        primary = np.ones_like(gfracmasked, dtype='?')
    mws = primary.copy()

    # ADM apply the mask/logic selection for all MWS-MAIN targets
    mws &= (gfracmasked < 0.5) &(gflux > 0) & (gnobs > 0)
    mws &= (rfracmasked < 0.5) & (rflux > 0) & (rnobs > 0)

    mws &= (objtype == b'PSF')
    return mws

############################################################
def isMWS_main(gflux=None, rflux=None, zflux=None, w1flux=None, w2flux=None,
               gnobs=None, rnobs=None, gfracmasked=None, rfracmasked=None,
               pmra=None, pmdec=None, parallax=None, obs_rflux=None,
               gaia=None, gaiagmag=None, gaiabmag=None, gaiarmag=None,
               gaiaaen=None, gaiadupsource=None, 
               objtype=None, primary=None, south=True):
    """Set bits for main ``MWS`` targets.

    Args:
        see :func:`~desitarget.cuts.set_target_bits` for parameters.

    Returns:
        mask1 : array_like.
            ``True`` if and only if the object is a ``MWS_BROAD`` target.
        mask2 : array_like.
            ``True`` if and only if the object is a ``MWS_MAIN_RED`` target.
        mask3 : array_like.
            ``True`` if and only if the object is a ``MWS_MAIN_BLUE`` target.

    Notes:
        - as of 11/2/18, based on version 158 on `the wiki`_.
    """
    from desitarget.cuts import notinMWS_main_mask
    from desitarget.cuts import isMWS_main_colors
    from desitarget.cuts import notinMWS_faint_mask
    from desitarget.cuts import isMWS_faint_colors
    
    if primary is None:
        primary = np.ones_like(gaia, dtype='?')
    mws = primary.copy()

    # ADM currently no difference between N/S for MWS, so easiest
    # ADM just to use one selection
    # if south:

    # ADM do not target any objects for which entries are NaN
    # ADM and turn off the NaNs for those entries
    nans = (np.isnan(rflux) | np.isnan(gflux) |
            np.isnan(parallax) | np.isnan(pmra) | np.isnan(pmdec))
    w = np.where(nans)[0]
    if len(w) > 0:
        # ADM make copies as we are reassigning values
        rflux, gflux, obs_rflux = rflux.copy(), gflux.copy(), obs_rflux.copy()
        parallax, pmra, pmdec = parallax.copy(), pmra.copy(), pmdec.copy()
        rflux[w], gflux[w], obs_rflux[w] = 0., 0., 0.
        parallax[w], pmra[w], pmdec[w] = 0., 0., 0.
        mws &= ~nans
        log.info('{}/{} NaNs in file...t = {:.1f}s'
                 .format(len(w), len(mws), time()-start))

    mws_main = notinMWS_main_mask(gaia=gaia, gfracmasked=gfracmasked, gnobs=gnobs,
                              gflux=gflux, rfracmasked=rfracmasked, rnobs=rnobs,
                              rflux=rflux, gaiadupsource=gaiadupsource, 
                              objtype=objtype, primary=primary)

    mws_faint = notinMWS_faint_mask(gfracmasked=gfracmasked, gnobs=gnobs,
                              gflux=gflux, rfracmasked=rfracmasked, rnobs=rnobs,
                              rflux=rflux, 
                              objtype=objtype, primary=primary)

    
    # ADM pass the mws that pass cuts as primary, to restrict to the
    # ADM sources that weren't in a mask/logic cut.
    mws_broad, mws_red, mws_blue = isMWS_main_colors(
        gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, w2flux=w2flux,
        pmra=pmra, pmdec=pmdec, parallax=parallax, obs_rflux=obs_rflux,
        gaiagmag=gaiagmag, gaiabmag=gaiabmag, gaiarmag=gaiarmag, gaiaaen=gaiaaen,
        primary=mws_main, south=south
    )

    mws_red_faint, mws_blue_faint = isMWS_faint_colors(
        gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, w2flux=w2flux,
        obs_rflux=obs_rflux,
        primary=mws_faint, south=south
    )
    
    return mws_broad, mws_red, mws_blue, mws_red_faint, mws_blue_faint

############################################################
def set_target_bits(photsys_north, photsys_south, obs_rflux,
                    gflux, rflux, zflux, w1flux, w2flux,
                    objtype, release, gfluxivar, rfluxivar, zfluxivar,
                    gnobs, rnobs, znobs, gfracflux, rfracflux, zfracflux,
                    gfracmasked, rfracmasked, zfracmasked,
                    gfracin, rfracin, zfracin, gallmask, rallmask, zallmask,
                    gsnr, rsnr, zsnr, w1snr, w2snr, deltaChi2, dchisq,
                    gaia, pmra, pmdec, parallax, parallaxovererror, parallaxerr,
                    gaiagmag, gaiabmag, gaiarmag, gaiaaen, gaiadupsource,
                    gaiaparamssolved, gaiabprpfactor, gaiasigma5dmax, galb,
                    tcnames, qso_optical_cuts, qso_selection, brightstarinblob,
                    Grr, primary):
    """Perform target selection on parameters, returning target mask arrays.

    Parameters
    ----------
    photsys_north, photsys_south : :class:`~numpy.ndarray`
        ``True`` for objects that were drawn from northern (MzLS/BASS) or
        southern (DECaLS) imaging, respectively.
    obs_rflux : :class:`~numpy.ndarray`
        `rflux` but WITHOUT any Galactic extinction correction.
    gflux, rflux, zflux, w1flux, w2flux : :class:`~numpy.ndarray`
        The flux in nano-maggies of g, r, z, W1 and W2 bands.
    objtype, release : :class:`~numpy.ndarray`
        `The Legacy Surveys`_ imaging ``TYPE`` and ``RELEASE`` columns.
    gfluxivar, rfluxivar, zfluxivar: :class:`~numpy.ndarray`
        The flux inverse variances in g, r, and z bands.
    gnobs, rnobs, znobs: :class:`~numpy.ndarray`
        The number of observations (in the central pixel) in g, r and z.
    gfracflux, rfracflux, zfracflux: :class:`~numpy.ndarray`
        Profile-weighted fraction of the flux from other sources divided
        by the total flux in g, r and z bands.
    gfracmasked, rfracmasked, zfracmasked: :class:`~numpy.ndarray`
        Fraction of masked pixels in the g, r and z bands.
    gallmask, rallmask, zallmask: :class:`~numpy.ndarray`
        Bitwise mask set if the central pixel from all images
        satisfy each condition in g, r, z.
    gsnr, rsnr, zsnr, w1snr, w2snr: :class:`~numpy.ndarray`
        Signal-to-noise in g, r, z, W1 and W2 defined as the flux per
        band divided by sigma (flux x sqrt of the inverse variance).
    deltaChi2: :class:`~numpy.ndarray`
        chi2 difference between PSF and SIMP, dchisq_PSF - dchisq_SIMP.
    dchisq: :class:`~numpy.ndarray`
        Difference in chi2  between successively more-complex model fits.
        Columns are model fits, in order, of PSF, REX, EXP, DEV, COMP.
    gaia: :class:`~numpy.ndarray`
        ``True`` if there is a match between this object in
        `the Legacy Surveys`_ and in Gaia.
    pmra, pmdec, parallax, parallaxovererror: :class:`~numpy.ndarray`
        Gaia-based proper motion in RA and Dec, and parallax and error.
    gaiagmag, gaiabmag, gaiarmag: :class:`~numpy.ndarray`
            Gaia-based g-, b- and r-band MAGNITUDES.
    gaiaaen, gaiadupsource, gaiaparamssolved: :class:`~numpy.ndarray`
        Gaia-based measures of Astrometric Excess Noise, whether the source
        is a duplicate, and how many parameters were solved for.
    gaiabprpfactor, gaiasigma5dmax: :class:`~numpy.ndarray`
        Gaia_based BP/RP excess factor and longest semi-major axis
        of 5-d error ellipsoid.
    galb: :class:`~numpy.ndarray`
        Galactic latitude (degrees).
    tcnames : :class:`list`, defaults to running all target classes
        A list of strings, e.g. ['QSO','LRG']. If passed, process targeting only
        for those specific target classes. A useful speed-up when testing.
        Options include ["ELG", "QSO", "LRG", "MWS", "BGS", "STD"].
    qso_optical_cuts : :class:`boolean` defaults to ``False``
        Apply just optical color-cuts when selecting QSOs with
        ``qso_selection="colorcuts"``.
    qso_selection : :class:`str`, optional, defaults to ``'randomforest'``
        The algorithm to use for QSO selection; valid options are
        ``'colorcuts'`` and ``'randomforest'``
    brightstarinblob: boolean array_like or None
        ``True`` if the object shares a blob with a "bright" (Tycho-2) star.
    Grr: array_like or None
        Gaia G band magnitude minus observational r magnitude.
    primary : :class:`~numpy.ndarray`
        ``True`` for objects that should be considered when setting bits.
    survey : :class:`str`, defaults to ``'main'``
        Specifies which target masks yaml file and target selection cuts
        to use. Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.

    Returns
    -------
    :class:`~numpy.ndarray`
        (desi_target, bgs_target, mws_target) where each element is
        an ndarray of target selection bitmask flags for each object.

    Notes
    -----
    - Gaia quantities have units that are the same as `the Gaia data model`_.
    """
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask
    from desitarget.cuts import isBGS
    from desitarget.cuts import isMWS_main, isMWS_nearby, isMWS_WD
    from desitarget.cuts import isELG, isLRGpass, isQSO_cuts, isQSO_randomforest
    from desitarget.cuts import isSTD
    
    if "LRG" in tcnames:
        lrg_north, lrg1pass_north, lrg2pass_north = isLRGpass(
            primary=primary,
            gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, gflux_ivar=gfluxivar,
            rflux_snr=rsnr, zflux_snr=zsnr, w1flux_snr=w1snr, south=False
        )

        lrg_south, lrg1pass_south, lrg2pass_south = isLRGpass(
            primary=primary,
            gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, gflux_ivar=gfluxivar,
            rflux_snr=rsnr, zflux_snr=zsnr, w1flux_snr=w1snr, south=True
        )
    else:
        # ADM if not running the LRG selection, set everything to arrays of False
        lrg_north, lrg1pass_north, lrg2pass_north = ~primary, ~primary, ~primary
        lrg_south, lrg1pass_south, lrg2pass_south = ~primary, ~primary, ~primary

    # ADM combine LRG target bits for an LRG target based on any imaging
    lrg = (lrg_north & photsys_north) | (lrg_south & photsys_south)
    lrg1pass = (lrg1pass_north & photsys_north) | (lrg1pass_south & photsys_south)
    lrg2pass = (lrg2pass_north & photsys_north) | (lrg2pass_south & photsys_south)

    if "ELG" in tcnames:
        elg_classes = []
        for south in [False, True]:
            elg_classes.append(
                isELG(
                    primary=primary, gflux=gflux, rflux=rflux, zflux=zflux,
                    gallmask=gallmask, rallmask=rallmask, zallmask=zallmask,
                    brightstarinblob=brightstarinblob, south=south)
            )
        elg_north, elg_south = elg_classes
    else:
        # ADM if not running the ELG selection, set everything to arrays of False.
        elg_north, elg_south = ~primary, ~primary

    # ADM combine ELG target bits for an ELG target based on any imaging
    elg = (elg_north & photsys_north) | (elg_south & photsys_south)

    if "QSO" in tcnames:
        if qso_selection == 'colorcuts':
            # ADM determine quasar targets in the north and the south separately
            qso_north = isQSO_cuts(
                primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
                w1flux=w1flux, w2flux=w2flux,
                deltaChi2=deltaChi2, brightstarinblob=brightstarinblob,
                objtype=objtype, w1snr=w1snr, w2snr=w2snr, release=release,
                optical=qso_optical_cuts, south=False
            )
            qso_south = isQSO_cuts(
                primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
                w1flux=w1flux, w2flux=w2flux,
                deltaChi2=deltaChi2, brightstarinblob=brightstarinblob,
                objtype=objtype, w1snr=w1snr, w2snr=w2snr, release=release,
                optical=qso_optical_cuts, south=True
            )
        elif qso_selection == 'randomforest':
            # ADM determine quasar targets in the north and the south separately
            qso_north = isQSO_randomforest(
                primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
                w1flux=w1flux, w2flux=w2flux,
                deltaChi2=deltaChi2, brightstarinblob=brightstarinblob,
                objtype=objtype, release=release, south=False
            )
            qso_south = isQSO_randomforest(
                primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
                w1flux=w1flux, w2flux=w2flux,
                deltaChi2=deltaChi2, brightstarinblob=brightstarinblob,
                objtype=objtype, release=release, south=True
            )
        else:
            raise ValueError('Unknown qso_selection {}; valid options are {}'.format(
                qso_selection, qso_selection_options))
    else:
        # ADM if not running the QSO selection, set everything to arrays of False
        qso_north, qso_south = ~primary, ~primary

    # ADM combine quasar target bits for a quasar target based on any imaging
    qso = (qso_north & photsys_north) | (qso_south & photsys_south)

    # ADM set the BGS bits
    if "BGS" in tcnames:
        bgs_classes = []
        for targtype in ["bright", "faint", "wise"]:
            for south in [False, True]:
                bgs_classes.append(
                    isBGS(
                        gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, w2flux=w2flux,
                        gnobs=gnobs, rnobs=rnobs, znobs=znobs,
                        gfracmasked=gfracmasked, rfracmasked=rfracmasked, zfracmasked=zfracmasked,
                        gfracflux=gfracflux, rfracflux=rfracflux, zfracflux=zfracflux,
                        gfracin=gfracin, rfracin=rfracin, zfracin=zfracin,
                        gfluxivar=gfluxivar, rfluxivar=rfluxivar, zfluxivar=zfluxivar,
                        brightstarinblob=brightstarinblob, Grr=Grr, w1snr=w1snr, gaiagmag=gaiagmag,
                        objtype=objtype, primary=primary, south=south, targtype=targtype
                    )
                )

        bgs_bright_north, bgs_bright_south,      \
            bgs_faint_north, bgs_faint_south,    \
            bgs_wise_north, bgs_wise_south =     \
                                                 bgs_classes
    else:
        # ADM if not running the BGS selection, set everything to arrays of False
        bgs_bright_north, bgs_bright_south = ~primary, ~primary
        bgs_faint_north, bgs_faint_south = ~primary, ~primary
        bgs_wise_north, bgs_wise_south = ~primary, ~primary

    # ADM combine BGS targeting bits for a BGS selected in any imaging
    bgs_bright = (bgs_bright_north & photsys_north) | (bgs_bright_south & photsys_south)
    bgs_faint = (bgs_faint_north & photsys_north) | (bgs_faint_south & photsys_south)
    bgs_wise = (bgs_wise_north & photsys_north) | (bgs_wise_south & photsys_south)

    if "MWS" in tcnames:
        mws_classes = []
        # ADM run the MWS target types for both north and south
        for south in [False, True]:
            mws_classes.append(
                isMWS_main(
                    gaia=gaia, gaiaaen=gaiaaen, gaiadupsource=gaiadupsource,
                    gflux=gflux, rflux=rflux, obs_rflux=obs_rflux,
                    gnobs=gnobs, rnobs=rnobs,
                    gfracmasked=gfracmasked, rfracmasked=rfracmasked,
                    pmra=pmra, pmdec=pmdec, parallax=parallax,
                    objtype=objtype,
                    primary=primary, south=south
                )
            )

        mws_broad_n, mws_red_n, mws_blue_n,  \
        mws_red_faint_n, mws_blue_faint_n,   \
        mws_broad_s, mws_red_s, mws_blue_s,  \
        mws_red_faint_s, mws_blue_faint_s =  np.vstack(mws_classes)

        mws_nearby = isMWS_nearby(
            gaia=gaia, gaiagmag=gaiagmag, parallax=parallax,
            parallaxerr=parallaxerr
        )
        
        mws_wd = isMWS_WD(
            gaia=gaia, galb=galb, astrometricexcessnoise=gaiaaen,
            pmra=pmra, pmdec=pmdec, parallax=parallax, parallaxovererror=parallaxovererror,
            photbprpexcessfactor=gaiabprpfactor, astrometricsigma5dmax=gaiasigma5dmax,
            gaiagmag=gaiagmag, gaiabmag=gaiabmag, gaiarmag=gaiarmag
        )
    else:
        # ADM if not running the MWS selection, set everything to arrays of False
        mws_broad_n, mws_red_n, mws_blue_n = ~primary, ~primary, ~primary
        mws_broad_s, mws_red_s, mws_blue_s = ~primary, ~primary, ~primary
        mws_nearby, mws_wd = ~primary, ~primary
        
        mws_blue_faint   = ~primary
        mws_blue_faint_s = ~primary
        mws_blue_faint_n = ~primary
        
        mws_red_faint   = ~primary
        mws_red_faint_s = ~primary
        mws_red_faint_n = ~primary

    if "STD" in tcnames:
        # ADM Make sure to pass all of the needed columns! At one point we stopped
        # ADM passing objtype, which meant no standards were being returned.
        std_faint = isSTD(
            primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
            gfracflux=gfracflux, rfracflux=rfracflux, zfracflux=zfracflux,
            gfracmasked=gfracmasked, rfracmasked=rfracmasked, objtype=objtype,
            zfracmasked=zfracmasked, gnobs=gnobs, rnobs=rnobs, znobs=znobs,
            gfluxivar=gfluxivar, rfluxivar=rfluxivar, zfluxivar=zfluxivar,
            gaia=gaia, astrometricexcessnoise=gaiaaen, paramssolved=gaiaparamssolved,
            pmra=pmra, pmdec=pmdec, parallax=parallax, dupsource=gaiadupsource,
            gaiagmag=gaiagmag, gaiabmag=gaiabmag, gaiarmag=gaiarmag, bright=False
        )
        std_bright = isSTD(
            primary=primary, zflux=zflux, rflux=rflux, gflux=gflux,
            gfracflux=gfracflux, rfracflux=rfracflux, zfracflux=zfracflux,
            gfracmasked=gfracmasked, rfracmasked=rfracmasked, objtype=objtype,
            zfracmasked=zfracmasked, gnobs=gnobs, rnobs=rnobs, znobs=znobs,
            gfluxivar=gfluxivar, rfluxivar=rfluxivar, zfluxivar=zfluxivar,
            gaia=gaia, astrometricexcessnoise=gaiaaen, paramssolved=gaiaparamssolved,
            pmra=pmra, pmdec=pmdec, parallax=parallax, dupsource=gaiadupsource,
            gaiagmag=gaiagmag, gaiabmag=gaiabmag, gaiarmag=gaiarmag, bright=True
        )
        # ADM the standard WDs are currently identical to the MWS WDs
        std_wd = mws_wd
    else:
        # ADM if not running the standards selection, set everything to arrays of False
        std_faint, std_bright, std_wd = ~primary, ~primary, ~primary

    # ADM combine the north/south MWS bits.
    mws_broad = (mws_broad_n & photsys_north) | (mws_broad_s & photsys_south)
    mws_blue  = (mws_blue_n & photsys_north) | (mws_blue_s & photsys_south)
    mws_red   = (mws_red_n & photsys_north) | (mws_red_s & photsys_south)

    # APC
    mws_blue_faint = (mws_blue_faint_n & photsys_north) | (mws_blue_faint_s & photsys_south)
    mws_red_faint  = (mws_red_faint_n  & photsys_north) | (mws_red_faint_s  & photsys_south)
    
    # Construct the targetflag bits for DECaLS (i.e. South).
    desi_target = lrg_south * desi_mask.LRG_SOUTH
    desi_target |= elg_south * desi_mask.ELG_SOUTH
    desi_target |= qso_south * desi_mask.QSO_SOUTH

    # Construct the targetflag bits for MzLS and BASS (i.e. North).
    desi_target |= lrg_north * desi_mask.LRG_NORTH
    desi_target |= elg_north * desi_mask.ELG_NORTH
    desi_target |= qso_north * desi_mask.QSO_NORTH

    # Construct the targetflag bits combining north and south.
    desi_target |= lrg * desi_mask.LRG
    desi_target |= elg * desi_mask.ELG
    desi_target |= qso * desi_mask.QSO

    # ADM add the per-pass information in the south...
    desi_target |= lrg1pass_south * desi_mask.LRG_1PASS_SOUTH
    desi_target |= lrg2pass_south * desi_mask.LRG_2PASS_SOUTH
    # ADM ...the north...
    desi_target |= lrg1pass_north * desi_mask.LRG_1PASS_NORTH
    desi_target |= lrg2pass_north * desi_mask.LRG_2PASS_NORTH
    # ADM ...and combined.
    desi_target |= lrg1pass * desi_mask.LRG_1PASS
    desi_target |= lrg2pass * desi_mask.LRG_2PASS

    # ADM Standards.
    desi_target |= std_faint * desi_mask.STD_FAINT
    desi_target |= std_bright * desi_mask.STD_BRIGHT
    desi_target |= std_wd * desi_mask.STD_WD

    # BGS bright and faint, south.
    bgs_target = bgs_bright_south * bgs_mask.BGS_BRIGHT_SOUTH
    bgs_target |= bgs_faint_south * bgs_mask.BGS_FAINT_SOUTH
    bgs_target |= bgs_wise_south * bgs_mask.BGS_WISE_SOUTH

    # BGS bright and faint, north.
    bgs_target |= bgs_bright_north * bgs_mask.BGS_BRIGHT_NORTH
    bgs_target |= bgs_faint_north * bgs_mask.BGS_FAINT_NORTH
    bgs_target |= bgs_wise_north * bgs_mask.BGS_WISE_NORTH

    # BGS combined, bright and faint
    bgs_target |= bgs_bright * bgs_mask.BGS_BRIGHT
    bgs_target |= bgs_faint * bgs_mask.BGS_FAINT
    bgs_target |= bgs_wise * bgs_mask.BGS_WISE

    # ADM MWS main, nearby, and WD.
    mws_target = mws_broad * mws_mask.MWS_BROAD
    mws_target |= mws_wd * mws_mask.MWS_WD
    mws_target |= mws_nearby * mws_mask.MWS_NEARBY

    # ADM MWS main north/south split.
    mws_target |= mws_broad_n * mws_mask.MWS_BROAD_NORTH
    mws_target |= mws_broad_s * mws_mask.MWS_BROAD_SOUTH
    
    # ADM MWS main blue/red split.
    mws_target |= mws_blue * mws_mask.MWS_MAIN_BLUE
    mws_target |= mws_blue_n * mws_mask.MWS_MAIN_BLUE_NORTH
    mws_target |= mws_blue_s * mws_mask.MWS_MAIN_BLUE_SOUTH
    mws_target |= mws_red * mws_mask.MWS_MAIN_RED
    mws_target |= mws_red_n * mws_mask.MWS_MAIN_RED_NORTH
    mws_target |= mws_red_s * mws_mask.MWS_MAIN_RED_SOUTH

    # APC
    mws_target |= mws_blue_faint   * mws_mask.MWS_MAIN_BLUE_FAINT
    mws_target |= mws_blue_faint_n * mws_mask.MWS_MAIN_BLUE_FAINT_NORTH
    mws_target |= mws_blue_faint_s * mws_mask.MWS_MAIN_BLUE_FAINT_SOUTH

    mws_target |= mws_red_faint   * mws_mask.MWS_MAIN_RED_FAINT
    mws_target |= mws_red_faint_n * mws_mask.MWS_MAIN_RED_FAINT_NORTH
    mws_target |= mws_red_faint_s * mws_mask.MWS_MAIN_RED_FAINT_SOUTH

    # Are any BGS or MWS bit set?  Tell desi_target too.
    desi_target |= (bgs_target != 0) * desi_mask.BGS_ANY
    desi_target |= (mws_target != 0) * desi_mask.MWS_ANY

    return desi_target, bgs_target, mws_target