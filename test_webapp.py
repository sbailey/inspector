import os
import unittest
from io import BytesIO
import warnings
import base64
from astropy.table import Table
import numpy as np

from app import app

#- Utility to temporarily override environment variables while safely cleaning up
#- Written by LBL CBorg Coder AI
from contextlib import ContextDecorator
class TempEnvironmentVariable(ContextDecorator):
    def __init__(self, key, value):
        self.key = key
        self.value = value
        self.original_value = None

    def __enter__(self):
        self.original_value = os.getenv(self.key)
        if self.value is None:
            if self.original_value is not None:
                del os.environ[self.key]
        else:
            os.environ[self.key] = self.value

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.original_value is not None:
            os.environ[self.key] = self.original_value
        elif self.key in os.environ:
            del os.environ[self.key]


class FlaskAppTestCase(unittest.TestCase):
    def setUp(self):
        self.app = app.test_client()
        self.app.testing = True

        warnings.filterwarnings("ignore", message=r".*Cannot merge meta key .*", category=Warning)

    def tearDown(self):
        warnings.resetwarnings()

    #---------------------------------------------------------------------
    #- non-data pages

    def test_info_pages(self):
        for page in ('', 'search', 'examples', 'license', 'about', 'robots.txt', 'dr1/testauth'):
            response = self.app.get(f'/{page}')
            self.assertEqual(response.status_code, 200)

    def test_nonexistent_page(self):
        response = self.app.get('/nonexistent')
        self.assertEqual(response.status_code, 404)

    def test_auth(self):
        """Confirm that all proprietary endpoints require authorization"""
        #- these should give 401 = unauthorized
        for release in ('dr2', 'daily', 'blatfoo'):
            for targspec in ('targets', 'spectra'):
                for tilepix in ('', 'tiles', 'healpix'):
                    for search in ('radec/210,5,30',  # RA,dec cone search
                                   '150/0-5',         # TILEID/FIBERs
                                   '1234,5678'        # TARGETIDs
                                   ):
                        #- TILEID/FIBERs search doesn't apply for healpix
                        if tilepix=='healpix' and search=='150/0-5':
                            continue

                        url = f'/{release}/{targspec}/{tilepix}/{search}'.replace('//', '/')
                        response = self.app.get(url)
                        errmsg = f"{url} didn't return 401=Unauthorized"
                        self.assertEqual(response.status_code, 401, errmsg)

        #- successful authorization with fake credentials
        with TempEnvironmentVariable('DESI_COLLAB_USERNAME', 'blat'):
            with TempEnvironmentVariable('DESI_COLLAB_PASSWORD', 'foo'):
                #- success
                auth_bytes = "blat:foo".encode('utf-8')
                base64_bytes = base64.b64encode(auth_bytes)
                base64_string = base64_bytes.decode('utf-8')
                headers = {'Authorization': f'Basic {base64_string}'}
                response = self.app.get('/dr2/targets/tiles/7614/1500-1503', headers=headers)
                self.assertEqual(response.status_code, 200)

                #- wrong password = fail
                auth_bytes = "blat:bizbat".encode('utf-8')
                base64_bytes = base64.b64encode(auth_bytes)
                base64_string = base64_bytes.decode('utf-8')
                headers = {'Authorization': f'Basic {base64_string}'}
                response = self.app.get('/dr2/targets/tiles/7614/1500-1503', headers=headers)
                self.assertEqual(response.status_code, 401)



    #---------------------------------------------------------------------
    #- Targets

    def test_targets_radec(self):
        #- basic
        response = self.app.get('/dr1/targets/radec/210,5,30')
        self.assertEqual(response.status_code, 200)

        #- healpix is default
        response_healpix = self.app.get('/dr1/targets/healpix/radec/210,5,30')
        self.assertEqual(response_healpix.status_code, 200)
        data = response_healpix.data.replace(b'/healpix/', b'/')
        self.assertEqual(data, response.data)

        #- tiles
        response_tiles = self.app.get('/dr1/targets/tiles/radec/210,5,30')
        self.assertEqual(response_tiles.status_code, 200)
        self.assertNotEqual(response_healpix.data, response_tiles.data)

        #- radius optional if within 10" of a target
        response = self.app.get('/dr1/targets/radec/210,4.9946')
        self.assertEqual(response.status_code, 200)

        #- other data releases
        for release in ('edr', 'fuji', 'iron'):
            response = self.app.get(f'/{release}/targets/radec/210,4.9946')
            self.assertEqual(response.status_code, 200)

    def test_targets_radec_formats(self):
        for format in ('html', 'fits', 'json'):
            response = self.app.get(f'/dr1/targets/radec/210,5,30?format={format}')
            self.assertEqual(response.status_code, 200)

        #- csv and ascii should be directly parseable by astropy.table.Table
        for format in ('csv', 'ascii'):
            response = self.app.get(f'/dr1/targets/radec/210,5,30?format={format}')
            t = Table.read(BytesIO(response.data), format=format)
            self.assertEqual(len(t), 4)

    def test_targets_radec_failures(self):
        #- 404 Not Found if no targets found, e.g. outside of footprint
        response = self.app.get('/dr1/targets/radec/210,-80,10')
        self.assertEqual(response.status_code, 404)

        #- 400 Bad Request if ra,dec out of range
        response = self.app.get('/dr1/targets/radec/210,100,10')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/targets/radec/210,-100,10')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/targets/radec/410,10,10')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/targets/radec/-10,10,10')
        self.assertEqual(response.status_code, 400)

        #- malformed ra,dec,radius
        response = self.app.get('/dr1/targets/radec/210')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/targets/radec/1,2,3,4')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/targets/radec/200deg,10deg,5arcsec')
        self.assertEqual(response.status_code, 400)

        #- invalid format
        response = self.app.get('/dr1/targets/radec/210,5,30?format=invalid')
        self.assertEqual(response.status_code, 400)

    def test_targets_tilefibers(self):
        #- basic
        response = self.app.get('/dr1/targets/150/0-5,10,20:22')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/targets/tiles/150/0-5,10,20:22')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/targets/tiles/150/0-5,10,20:22?format=fits')
        self.assertEqual(response.status_code, 200)

        #- various ways of specifying fiber ranges
        t = Table.read(BytesIO(self.app.get('/dr1/targets/tiles/150/0:5?format=csv').data), format='csv')
        self.assertEqual(len(t), 5)
        self.assertTrue(np.all(t['FIBER'] == np.arange(5)))

        t = Table.read(BytesIO(self.app.get('/dr1/targets/tiles/150/0-5?format=csv').data), format='csv')
        self.assertEqual(len(t), 6)
        self.assertTrue(np.all(t['FIBER'] == np.arange(6)))

        t = Table.read(BytesIO(self.app.get('/dr1/targets/tiles/150/10,5,0?format=csv').data), format='csv')
        self.assertEqual(len(t), 3)
        self.assertTrue(np.all(t['FIBER'] == np.array([10,5,0])))

    def test_targets_tiles_failures(self):
        #- tile not found
        #- 404 Not Found; not 400 Bad Request, since someday that tile might exists
        response = self.app.get('/dr1/targets/tiles/999999/0-5')
        self.assertEqual(response.status_code, 404)

        #- fibers out of range (400 Bad Request because these fibers will never be valid)
        response = self.app.get('/dr1/targets/tiles/150/10000-10002')
        self.assertEqual(response.status_code, 400)

        #- bad format
        response = self.app.get('/dr1/targets/150/0-5?format=invalid')
        self.assertEqual(response.status_code, 400)

    def test_targets_targetids(self):
        #- basic
        response = self.app.get('/dr1/targets/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/targets/healpix/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/targets/tiles/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)

        t = Table.read(BytesIO(self.app.get('/dr1/targets/39627908959964170,39627908959964322?format=csv').data), format='csv')
        self.assertEqual(len(t), 3)  #- one target observed in both sv3 and main

        #- filter by SURVEY
        t = Table.read(BytesIO(self.app.get('/dr1/targets/39627908959964170,39627908959964322?SURVEY=sv3&format=csv').data), format='csv')
        self.assertEqual(len(t), 2)


    #---------------------------------------------------------------------
    #- Spectra
    #- slower, so less complete coverage of all matrix options

    def test_spectra_radec(self):
        #- basic
        response = self.app.get('/dr1/spectra/radec/210,5,30')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/spectra/healpix/radec/210,5,30')
        self.assertEqual(response.status_code, 200)

        response = self.app.get('/dr1/spectra/tiles/radec/210,5,30')
        self.assertEqual(response.status_code, 200)

        #- fits
        response = self.app.get('/dr1/spectra/radec/210,5,30?format=fits')
        self.assertEqual(response.status_code, 200)

        #- plotting with noise model
        response = self.app.get('/dr1/spectra/radec/210,5,30?plotnoise=1')
        self.assertEqual(response.status_code, 200)

    def test_spectra_tilefibers(self):
        #- basic
        response = self.app.get('/dr1/spectra/150/0,2,3')
        self.assertEqual(response.status_code, 200)

    def test_spectra_targetids(self):
        #- basic (defaults to healpix)
        response = self.app.get('/dr1/spectra/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)
        #- healpix
        response = self.app.get('/dr1/spectra/healpix/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)
        #- tiles
        response = self.app.get('/dr1/spectra/tiles/39627908959964170,39627908959964322')
        self.assertEqual(response.status_code, 200)

    def test_spectra_failures(self):
        #- bad format
        response = self.app.get('/dr1/spectra/radec/210,5,30?format=invalid')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/spectra/tiles/1000/0-10?format=blatfoo')
        self.assertEqual(response.status_code, 400)

        #- no spectra found
        response = self.app.get('/dr1/spectra/radec/10,-80,10')
        self.assertEqual(response.status_code, 400)

        response = self.app.get('/dr1/spectra/1,2,3')
        self.assertEqual(response.status_code, 400)

        #- too many spectra for a single request
        response = self.app.get('/dr1/spectra/tiles/1000/0:5000')
        self.assertEqual(response.status_code, 400)

if __name__ == '__main__':
    unittest.main()
