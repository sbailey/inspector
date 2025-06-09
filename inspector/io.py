"""
inspector.io
============
"""

import numpy as np
from desispec import inventory
from desispec.io import read_spectra_parallel
from desispec.io.redrock import read_redrock_targetcat

MAX_RADIUS=1800  # 0.5 deg
MAX_RADIUS_ERROR_MESSAGE = f'Please limit your search to radius < {MAX_RADIUS} arcsec'

MAX_SPECTRA=1000
MAX_SPECTRA_ERROR_MESSAGE = '{} spectra is more than we can realistically display; please limit your search to fewer than {} spectra'

def standardize_specprod(specprod):
    """
    Return standardized specprod with aliases for dr1 -> iron, etc.
    """
    specprod_lower = specprod.lower()
    if specprod_lower == 'edr':
        return 'fuji'
    elif specprod_lower == 'dr1':
        return 'iron'
    elif specprod_lower == 'dr2':
        return 'loa'
    else:
        return specprod


def parse_fibers(fibers_string):
    """
    Parse fiber strings like "0-5,10:12,25" -> [0,1,2,3,4,5,10,11,25]
    Note: A-B is inclusive of B, while A:B is python-like exclusive of B
    """
    fibers = list()
    for token in fibers_string.split(','):
        if token.isdigit():
            fibers.append(int(token))
        elif '-' in token:
            first, last = map(int, token.split('-'))
            if first > last:
                raise ValueError(f'Failed to parse {token} as part of {fibers_string}')
            fibers.extend( range(first, last+1) )
        elif ':' in token:
            first, last = map(int, token.split(':'))
            if first > last:
                raise ValueError(f'Failed to parse {token} as part of {fibers_string}')
            fibers.extend( range(first, last) )
        else:
            raise ValueError(f'Failed to parse {token} as part of {fibers_string}')

    return fibers


def validate_radec(radec):
    ra,dec,radius = inventory.parse_radec(radec)
    if radius > MAX_RADIUS:
        raise ValueError(MAX_RADIUS_ERROR_MESSAGE)
    if ra<0 or ra>360:
        raise ValueError(f'RA={ra} should be between 0 and 360 degrees')
    if dec<-90 or dec>90:
        raise ValueError(f'dec={dec} should be between -90 and 90 degrees')

    return ra,dec,radius

def add_zcat_columns(targetcat, specprod):
    """
    Return copy of targetcat with zcat columns added
    
    Adds TARGET_RA, TARGET_DEC, SPECTYPE, Z, ZWARN
    """
    t = targetcat.copy(copy_data=False)
    specprod = standardize_specprod(specprod)
    zcat = read_redrock_targetcat(t, fmcols=['TARGET_RA', 'TARGET_DEC'], specprod=specprod)
    for col in ['TARGET_RA', 'TARGET_DEC', 'SPECTYPE', 'Z', 'ZWARN']:
        t[col] = zcat[col]

    if 'TARGETID' not in t.colnames:
        t['TARGETID'] = zcat['TARGETID']
    
    return t

def filter_table(table, filters):
    """
    Filter a Table based upon URL args dict

    TODO: expand explanation, args comes from request.args.to_dict(flat=False)
    """

    #- support blank or None filters
    if filters is None or len(filters) == 0:
        return table

    keep = np.ones(len(table), dtype=bool)
    for column, filter_list in filters.items():
        #- Interpret UPPERCASE as columns for filter; otherwise other options
        if not column.isupper():
            continue

        #- If UPPERCASE but not in the table columns, that is an error
        if column not in table.colnames:
            raise ValueError(f'Filter column "{column}" not in {table.colnames}')

        filter_list = np.atleast_1d(filter_list)
        for filt in filter_list:
            if isinstance(filt, str) and ':' in filt:
                operator, value = filt.split(':')
            else:
                operator = 'eq'
                value = filt

            value = np.array(value, dtype=table[column].dtype)
            ### print(f'Filter {column} {operator} {value}')

            if operator == 'eq':
                keep &= (table[column] == value)
            elif operator == 'lt':
                keep &= (table[column] < value)
            elif operator == 'le':
                keep &= (table[column] <= value)
            elif operator == 'gt':
                keep &= (table[column] > value)
            elif operator == 'ge':
                keep &= (table[column] >= value)
            elif operator == 'ne':
                keep &= (table[column] != value)
            else:
                #- TODO: handle error message with error message to user
                raise ValueError(f'Unrecognized operator {operator}')

    return table[keep]

def load_targets(specprod, specgroup, radec=None, targetids=None, filters=None):
    """
    required: specprod, specgroup; plus radec OR targetids (but not both)
    """
    specprod = standardize_specprod(specprod)

    if radec is not None:
        radec = validate_radec(radec)
    elif targetids is not None:
        if isinstance(targetids, str):
            targetids = list(map(int, targetids.split(',')))
    else:
        raise ValueError('must specify radec or targetids')

    if specgroup == 'healpix':
        t = inventory.target_healpix(radec=radec, targetids=targetids, specprod=specprod)
    elif specgroup == 'tiles':
        t = inventory.target_tiles(radec=radec, targetids=targetids, specprod=specprod)

    if len(t) > 0:
        t = add_zcat_columns(t, specprod)
        t = filter_table(t, filters)

    return t

def load_spectra(specprod, specgroup, radec=None, targetids=None, filters=None, maxspectra=MAX_SPECTRA):
    """
    Required: specprod, specgroup; and radec OR targetids (but not both)

    TODO: separate loading spectra from flask-specific rendering
    """
    specprod = standardize_specprod(specprod)
    targetcat = load_targets(specprod, specgroup=specgroup, radec=radec, targetids=targetids, filters=filters)

    num_spectra = len(targetcat)
    if num_spectra == 0:
        return None
    elif num_spectra > maxspectra:
        raise ValueError(MAX_SPECTRA_ERROR_MESSAGE.format(num_spectra, maxspectra))

    print(f'Reading {num_spectra} spectra')
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))
    return spectra

