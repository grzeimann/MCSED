""" MCSED - config.py


1) Configuration Settings

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
# Fit for metallicity
fit_metallicity = False
# Fixed metallicity of SSP models if fit_metallicity is False
metallicity = 0.0006  # for fixed metallicity

# SSP code for models
ssp = 'fsps'  # 'fsps'
isochrone = 'padova'  # 'padova', 'basti', 'mist', 'geneva', 'parsec'
sfh = 'double_powerlaw'
dust_law = 'noll'

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
             46: 'iracch4.res', 47: 'SubIB738.res', 48: 'SubIB797.res'}

# Catalog column name of filter and dictionary value to the filter file
catalog_filter_dict = {}
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

# Filters for testing code
test_filter_dict = {1: 'ia427', 2: 'ia445', 6: 'ia505', 7: 'ia527', 8: 'ia550',
                    9: 'ia574', 10: 'ia598', 11: 'ia624', 12: 'ia651',
                    13: 'ia679', 15: 'ia767', 17: 'ia856', 19: 'f435w',
                    20: 'f606wcand',  21: 'f814wcand', 22: 'f125w',
                    23: 'f140w', 24: 'f160w', 26: 'irac1', 27: 'irac2',
                    28: 'u', 29: 'r', 30: 'u38', 31: 'b', 32: 'v', 34: 'f606w',
                    35: 'rc', 36: 'f775w', 37: 'i', 38: 'f850lp',
                    39: 'f850lpcand', 40: 'j', 41: 'tenisj', 42: 'h',
                    43: 'tenisk', 44: 'ks', 45: 'irac3', 46: 'irac4',
                    47: 'ia738', 48: 'ia797'}
