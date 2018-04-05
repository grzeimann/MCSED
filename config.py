""" MCSED - config.py


1) Configuration Settings

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
# SSP code for models
ssp = 'fsps'  # options include: 'fsps'
isochrone = 'padova'  # options include: 'padova'
# SFH options include: 'double_powerlaw', 'empirical_direct', 'constant'
sfh = 'double_powerlaw'  
dust_law = 'noll'  # options include: 'noll', 'calzetti'
dust_em = 'DL07' # options include: 'DL07'

# Fix dust emission parameters
fix_dust_em = True

# EMCEE parameters
nwalkers = 100
nsteps = 1000

# Number of test objects
nobjects = 5

# Nebular Emission Properties
add_nebular = True
logU = -2.

# Error floor under which we don't trust the error estimates (fraction)
floor_error = 0.10

# Use metallicity-mass relationship from Ma et al. 2016
metallicity_mass_relationship = False
# Fixed metallicity of SSP models if fit_metallicity is False
fix_metallicity = False
metallicity = 0.0077  # for fixed metallicity

# Dictionaries
metallicity_dict = {'padova': [0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008,
                               0.0010, 0.0012, 0.0016, 0.0020, 0.0025, 0.0031,
                               0.0039, 0.0049, 0.0061, 0.0077, 0.0096, 0.0120,
                               0.0150, 0.0190, 0.0240, 0.0300]}

# FILTERS

# Common filter file names in FILTERS folder
filter_matrix_name = 'standard_filter_matrix.txt'
filt_dict = {0: 'SubB.res', 1: 'SubIB427.res', 2: 'SubIB445.res',
             3: 'SubIB464.res', 4: 'CFHTi.res', 5: 'SubIB484.res',
             6: 'SubIB505.res', 7: 'SubIB527.res', 8: 'SubIB550.res',
             9: 'SubIB574.res', 10: 'SubIB598.res', 11: 'SubIB624.res',
             12: 'SubIB651.res', 13: 'SubIB679.res', 14: 'SubIB709.res',
             15: 'SubIB767.res', 16: 'SubIB827.res',  17: 'SubIB856.res',
             18: 'SubIB907.res', 19: 'acs_f435w', 20: 'acs_f606w',
             21: 'acs_f814w', 22: 'wfc3_f125w', 23: 'wfc3_f140w',
             24: 'wfc3_f160w', 25: 'CFHTK.res', 26: 'iracch1.res',
             27: 'iracch2.res', 28: 'VLT_vimos_U.res', 29: 'VLT_vimos_R.res',
             30: 'lasilla22_WFI_U38.res', 31: 'lasilla22_WFI_B.res',
             32: 'lasilla22_WFI_V.res', 33: 'lasilla22_WFI_B.res',
             34: 'acs_f606w', 35: 'lasilla22_WFI_Rc.res', 36: 'acs_f775w',
             37: 'lasilla22_WFI_I.res', 38: 'acs_f850lp', 39: 'acs_f850lp',
             40: 'VLT_issac_J.res', 41: 'wircam_J.res', 42: 'VLT_issac_H.res',
             43: 'wircam_Ks.res', 44: 'VLT_issac_Ks.res', 45: 'iracch3.res',
             46: 'iracch4.res', 47: 'SubIB738.res', 48: 'SubIB797.res',
             49: 'KPNO_mosaicU.res', 50: 'keck_LRIS_g.res',
             51: 'keck_LRIS_Rs.res', 52: 'SubV.res', 53: 'Subrp.res',
             54: 'Subip.res', 55: 'Subzp.res', 56: 'Sub_moirics_H.res',
             57: 'Sub_moirics_J.res', 58: 'Sub_moirics_Ks.res',
             59: 'CFHTu.res', 60: 'CFHTg.res', 61: 'CFHTr.res',
             62: 'CFHTz.res', 63: 'ukidss_y', 64: 'ukidss_h', 65: 'ukidss_j',
             66: 'ukidss_k', 67: 'newfirm_J1.res', 68: 'newfirm_J2.res',
             69: 'newfirm_J3.res', 70: 'newfirm_H1.res', 71: 'newfirm_H2.res',
             72: 'newfirm_Ks.res', 73: 'wircam_H.res', 74: 'wircam_J.res',
             75: 'wircam_Ks.res', 76: 'MIPS24um.res'}

# Catalog column name of filter and dictionary value to the filter file
catalog_filter_dict, catalog_maglim_dict = {}, {}
catalog_filter_dict['goodss'] = {1: 'ia427', 2: 'ia445', 6: 'ia505',
                                 7: 'ia527', 8: 'ia550', 9: 'ia574',
                                 10: 'ia598', 11: 'ia624', 12: 'ia651',
                                 13: 'ia679', 15: 'ia767', 17: 'ia856',
                                 19: 'f435w', 20: 'f606wcand', 21: 'f814wcand',
                                 22: 'f125w', 23: 'f140w', 24: 'f160w',
                                 26: 'irac1', 27: 'irac2', 28: 'u', 29: 'r',
                                 30: 'u38', 31: 'b', 32: 'v', 34: 'f606w',
                                 35: 'rc', 36: 'f775w', 37: 'i', 38: 'f850lp',
                                 39: 'f850lpcand', 40: 'j', 41: 'tenisj',
                                 42: 'h', 43: 'tenisk', 44: 'ks', 45: 'irac3',
                                 46: 'irac4', 47: 'ia738', 48: 'ia797'}

catalog_filter_dict['goodsn'] = {0: 'b', 19: 'f435w', 22: 'f125w', 23: 'f140w',
                                 24: 'f160w', 26: 'irac1', 27: 'irac2',
                                 34: 'f606w', 36: 'f775w', 38: 'f850lp',
                                 45: 'irac3', 46: 'irac4', 49: 'u', 50: 'g',
                                 51: 'rs', 52: 'v', 53: 'r', 54: 'i', 55: 'z',
                                 56: 'h', 57: 'j', 58: 'ks'}

catalog_maglim_dict['goodsn'] = {0: 26.7, 19: 27.1, 22: 26.7, 23: 25.9,
                                 24: 26.1, 26: 24.5, 27: 24.6,
                                 34: 27.4, 36: 26.9, 38: 26.7,
                                 45: 22.8, 46: 22.7, 49: 26.4, 50: 26.3,
                                 51: 25.6, 52: 27.0, 53: 26.2, 54: 25.8,
                                 55: 25.5, 56: 24.3, 57: 25.0, 58: 24.7}

catalog_filter_dict['cosmos'] = {0: 'b', 1: 'ia427', 3: 'ia464', 4: 'i',
                                 5: 'ia484', 6: 'ia505', 7: 'ia527',
                                 9: 'ia574', 11: 'ia624', 13: 'ia679',
                                 14: 'ia709', 15: 'ia767', 16: 'ia827',
                                 20: 'f606w', 21: 'f814w', 22: 'f125w',
                                 23: 'f140w', 24: 'f160w', 26: 'irac1',
                                 27: 'irac2', 45: 'irac3', 46: 'irac4',
                                 47: 'ia738', 52: 'v', 53: 'rp', 54: 'ip',
                                 55: 'zp', 59: 'u', 60: 'g', 61: 'r', 62: 'z',
                                 63: 'uvista_y', 64: 'uvista_h',
                                 65: 'uvista_j', 66: 'uvista_ks', 67: 'j1',
                                 68: 'j2', 69: 'j3', 70: 'h1', 71: 'h2',
                                 72: 'k', 73: 'h', 74: 'j', 75: 'ks'}

catalog_maglim_dict['cosmos'] = {0: 27.9, 1: 26.3, 3: 25.8, 4: 27.0,
                                 5: 26.4, 6: 26.2, 7: 26.5,
                                 9: 25.7, 11: 26.4, 13: 25.9,
                                 14: 26.1, 15: 25.7, 16: 25.6,
                                 20: 26.7, 21: 26.5, 22: 26.1,
                                 23: 25.5, 24: 25.8, 26: 25.1,
                                 27: 25.0, 45: 21.6, 46: 21.6,
                                 47: 26.0, 52: 27.2, 53: 27.3, 54: 27.0,
                                 55: 25.8, 59: 26.7, 60: 27.2, 61: 27.2,
                                 62: 26.0,
                                 63: 25.7, 64: 24.9,
                                 65: 25.4, 66: 25.0, 67: 24.7,
                                 68: 24.5, 69: 24.4, 70: 23.7, 71: 23.6,
                                 72: 23.7, 73: 24.0, 74: 23.8, 75: 23.8}

catalog_filter_dict['aegis'] = {4: 'i', 20: 'f606w', 21: 'f814w', 22: 'f125w',
                                23: 'f140w', 24: 'f160w', 26: 'irac1',
                                27: 'irac2', 45: 'irac3', 46: 'irac4',
                                59: 'u', 60: 'g', 61: 'r', 62: 'z',
                                67: 'j1', 68: 'j2', 69: 'j3', 70: 'h1',
                                71: 'h2', 72: 'k', 73: 'h', 74: 'j', 75: 'ks'}

catalog_maglim_dict['aegis'] = {4: 27.0, 20: 26.8, 21: 26.4, 22: 26.3,
                                23: 25.7, 24: 26.1, 26: 25.2,
                                27: 25.1, 45: 22.6, 46: 22.6,
                                59: 26.8, 60: 27.4, 61: 27.3, 62: 26.2,
                                67: 24.8, 68: 24.5, 69: 24.4, 70: 23.8,
                                71: 23.8, 72: 23.6, 73: 24.3, 74: 24.5,
                                75: 24.0}
