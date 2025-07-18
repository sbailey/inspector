"""
Test some desispec.io functionality not covered by the webapp cases
"""

import unittest
import numpy as np
from astropy.table import Table

class TestIO(unittest.TestCase):

    def test_parse_fibers(self):
        from inspector.io import parse_fibers

        self.assertEqual(parse_fibers('10,20,30'), [10,20,30])
        self.assertEqual(parse_fibers('0:4'), [0,1,2,3])
        self.assertEqual(parse_fibers('0-4'), [0,1,2,3,4])
        self.assertEqual(parse_fibers('0:2,10,20-22'), [0,1,10,20,21,22])
        self.assertEqual(parse_fibers('100-103,0:2,10,20-22'), [100,101,102,103,0,1,10,20,21,22])

        with self.assertRaises(ValueError):
            fibers = parse_fibers('1,a,2')      # error: not integers

        with self.assertRaises(ValueError):
            fibers = parse_fibers('1,-3,2')     # error: negative fiber number

        with self.assertRaises(ValueError):
            fibers = parse_fibers('three-four') # error: not integers

        with self.assertRaises(ValueError):
            fibers = parse_fibers('10-5')       # error: first>last

        with self.assertRaises(ValueError):
            fibers = parse_fibers('10:5')       # error: first>last

    def test_validate_radec(self):
        from inspector.io import validate_radec

        self.assertEqual(validate_radec('10,20,30'), (10.0,20.0,30.0))
        self.assertEqual(validate_radec('10.5,-20,30'), (10.5,-20.0,30.0))
        self.assertEqual(validate_radec('10,-20.3'), (10.0,-20.3,10.0))

        with self.assertRaises(ValueError):
            ra,dec,radius = validate_radec('10,20,10000')    # radius too big

        with self.assertRaises(ValueError):
            ra,dec,radius = validate_radec('-10,20,30')      # RA<0

        with self.assertRaises(ValueError):
            ra,dec,radius = validate_radec('360.1,20,30')    # RA>360

        with self.assertRaises(ValueError):
            ra,dec,radius = validate_radec('10,-100,30')    # dec<-90

        with self.assertRaises(ValueError):
            ra,dec,radius = validate_radec('10,100,30')    # dec>90

    def test_filter_table(self):
        from astropy.table import Table
        from inspector.io import filter_table
        table = Table()
        table['A'] = list(range(10))

        #- scalar or list, number or string
        for value in [3, '3', [3,], ['3',], (3,)]:
            t = filter_table(table, dict(A=value))
            self.assertEqual(list(t['A']), [3,], f'Failed comparison test with {value=}')

        t = filter_table(table, dict(A='eq:5'))
        self.assertEqual(list(t['A']), [5,])

        t = filter_table(table, dict(A='gt:7'))
        self.assertEqual(list(t['A']), [8,9,])

        t = filter_table(table, dict(A='ge:7'))
        self.assertEqual(list(t['A']), [7,8,9,])

        t = filter_table(table, dict(A='lt:3'))
        self.assertEqual(list(t['A']), [0,1,2,])

        t = filter_table(table, dict(A='le:3'))
        self.assertEqual(list(t['A']), [0,1,2,3,])

        t = filter_table(table, dict(A='ne:4'))
        self.assertEqual(list(t['A']), [0,1,2,3,5,6,7,8,9])

        t = filter_table(table, dict(A=['ge:3', 'le:7']))
        self.assertEqual(list(t['A']), [3,4,5,6,7])

        with self.assertRaises(ValueError):
            t = filter_table(table, dict(A='xx:25'))

        with self.assertRaises(ValueError):
            t = filter_table(table, dict(A='eq:blat'))

    def test_load_targets(self):
        """Test load_targets, including failure modes"""
        from inspector.io import load_targets
        targets = load_targets('dr1', 'healpix', radec=(210,5,30))
        self.assertGreater(len(targets), 0)
        targets = load_targets('dr1', 'healpix', targetids=[39627908959964170, 39627908959964322])
        self.assertGreater(len(targets), 0)
        targets = load_targets('dr1', 'healpix', targetids='39627908959964170,39627908959964322')
        self.assertGreater(len(targets), 0)

        targets = load_targets('dr1', 'tiles', radec=(210,5,30))
        self.assertGreater(len(targets), 0)
        targets = load_targets('dr1', 'tiles', targetids=[39627908959964170, 39627908959964322])
        self.assertGreater(len(targets), 0)
        targets = load_targets('dr1', 'tiles', targetids='39627908959964170,39627908959964322')
        self.assertGreater(len(targets), 0)

        targets = load_targets('dr1', 'tiles', radec=(210,5,30), xcol=['SUBTYPE', 'FLUX_R'])
        self.assertIn('SUBTYPE', targets.colnames)
        self.assertIn('FLUX_R', targets.colnames)

        zerr_cut = 5e-5
        filters = dict(ZERR=f'gt:{zerr_cut}')
        targets = load_targets('dr1', 'tiles', radec=(210,5,30), filters=filters)
        self.assertIn('ZERR', targets.colnames)
        self.assertTrue(np.all(targets['ZERR']>zerr_cut))

        with self.assertRaises(ValueError):
            targets = load_targets('dr1', 'healpix')  # needs radec or targetids

    def test_load_spectra(self):
        """Test load_spectra (basic), including failure modes"""
        from inspector.io import load_spectra
        targetids = [39627908959964170, 39627908959964322]
        spectra = load_spectra('dr1', 'healpix', targetids=targetids)
        self.assertGreater(len(spectra), 0)

        spectra = load_spectra('dr1', 'healpix', targetids=[1,2,3])
        self.assertEqual(spectra, None)

        spectra = load_spectra('dr1', 'healpix', targetids=[])
        self.assertEqual(spectra, None)

        targetids = [39627908959964170, 39627908959964322]
        with self.assertRaises(ValueError):
            sp = load_spectra('dr1', 'healpix', targetids=targetids, maxspectra=1)

    def test_add_zcat_columns(self):
        from inspector.io import add_zcat_columns
        t1 = Table()
        t1 = Table()
        t1['FIBER'] = [0,1,2,3,4]
        t1['LASTNIGHT'] = 20210418
        t1['TILEID'] = 150

        t2 = add_zcat_columns(t1, 'iron')
        self.assertNotIn('Z', t1.colnames)  #- original wasn't changed
        self.assertIn('Z', t2.colnames)     #- new has extra columns
        self.assertIn('TARGETID', t2.colnames)

        #- New column in REDSHIFTS HDU
        t3 = add_zcat_columns(t1, 'iron', xcol=['SUBTYPE',])
        self.assertIn('SUBTYPE', t3.colnames)

        #- New columns in FIBERMAP HDU
        t4 = add_zcat_columns(t1, 'iron', xcol=['FLUX_G', 'FLUX_R'])
        self.assertIn('FLUX_G', t4.colnames)
        self.assertIn('FLUX_R', t4.colnames)

        #- xcol that was already in the list anyway
        t5 = add_zcat_columns(t1, 'iron', xcol=['Z', 'SPECTYPE'])
        self.assertEqual(t2.colnames, t5.colnames)

    def test_standardize_specprod(self):
        from inspector.io import standardize_specprod
        self.assertEqual(standardize_specprod('edr'), 'fuji')
        self.assertEqual(standardize_specprod('dr1'), 'iron')
        self.assertEqual(standardize_specprod('dr2'), 'loa')

        #- UPPERCASE DR names also work
        self.assertEqual(standardize_specprod('EDR'), 'fuji')
        self.assertEqual(standardize_specprod('DR1'), 'iron')
        self.assertEqual(standardize_specprod('DR2'), 'loa')

        #- Otherwise pass through specprod regardless of case
        self.assertEqual(standardize_specprod('blat'), 'blat')
        self.assertEqual(standardize_specprod('Blat'), 'Blat')


if __name__ == '__main__':
    unittest.main()
