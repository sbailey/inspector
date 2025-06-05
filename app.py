"""
DESI Data Inspector - app.py 
"""

import io
import os
import glob
import tempfile
from urllib.parse import urlencode
import hashlib
import uuid

import numpy as np
from astropy.table import Table
import fitsio

from desispec import inventory
from desispec.io import read_spectra_parallel, write_spectra, specprod_root
from desispec.io import findfile
from desispec.io.meta import get_lastnight
from desispec.io.redrock import read_redrock_targetcat

from prospect.viewer import plotspectra

from flask import Flask, request, jsonify, render_template, make_response, Response
import werkzeug.datastructures.structures

from inspector.auth import conditional_auth
from inspector.io import (standardize_specprod, parse_fibers, validate_radec,
                          load_targets, load_spectra,
                          filter_table, add_zcat_columns,
                          MAX_SPECTRA, MAX_SPECTRA_ERROR_MESSAGE)


app = Flask(__name__)
app.url_map.strict_slashes = False

### @app.route("/testargs")
### def test_filter():
###     print(type(request.args))
###     return str(request.args.to_dict(flat=False))

@app.route("/<string:specprod>/testauth")
@conditional_auth
def test_auth(specprod):
    return f"You are allowed to access {specprod} data"

@app.route("/")
def hello_world():
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    return render_template('homepage.html', root_url=root_url)

@app.route("/examples")
def examples():
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    return render_template('examples.html', root_url=root_url)

@app.route("/search")
def search():
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    return render_template('search.html', root_url=root_url)

@app.route("/license")
def license():
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    return render_template('license.html', root_url=root_url)

@app.route("/about")
def about():
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    return render_template('about.html', root_url=root_url)

@app.route("/robots.txt")
def robots():
    return Response("""User-agent: *
Disallow: /*
Allow: /
Allow: /examples
Allow: /search
Allow: /license
Allow: /about
""", mimetype="text/plain")

def _current_url_as_format(fmt):
    """
    Return current URL with alternate ?format=blah option
    """
    base_url = request.base_url  # the URL without the ?query=blat&format=foo options
    args = request.args.copy()
    args['format'] = fmt
    options = urlencode(args)
    newurl = f'{base_url}?{options}'
    return newurl

def render_table_html(table, header, description=''):
    """
    TODO: document
    """

    #- First check if table is empty
    if table is None or len(table) == 0:
        return render_template("error.html", code=404, summary='Not Found', message=description), 404

    #- And if there are too many targets
    ntargets = len(table)
    nmax = 10000
    if ntargets>nmax:
        msg = f"Your query resulted in {ntargets} targets, which is too many for us to realistically display."
        msg += f" Please limit your query to return less than {nmax} targets."
        return render_template("error.html", code=400, summary='Not Found', message=msg), 400

    #- limit precision for display
    for col in table.colnames:
        if table[col].dtype.kind == 'f':
            table[col].format = '{:.4f}'

    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    with io.StringIO() as buffer:
        table.write(buffer, format='ascii.html')
        table_html = buffer.getvalue()

    # construct download options footer text
    footer = 'Download table as'
    for fmt in ('json', 'csv', 'ascii', 'fits'):
        url = _current_url_as_format(fmt)
        footer += f' <a href={url}>{fmt}</a>'

    specview_url = _current_url_as_format('html').replace('/targets/', '/spectra/')
    specfits_url = _current_url_as_format('fits').replace('/targets/', '/spectra/')
    footer += f'; Spectra <a href="{specview_url}">view</a> <a href="{specfits_url}">fits</a>'

    if 'RA' in table.meta and 'DEC' in table.meta and 'RADIUS' in table.meta:
        ra = table.meta['RA']
        dec = table.meta['DEC']
        pixscale = table.meta['RADIUS']/128
    else:
        ra = dec = pixscale = None

    return render_template('inventory.html', root_url=root_url,
                           table_html=table_html,
                           header=header, description=description, footer=footer,
                           ra=ra, dec=dec, pixscale=pixscale)

def render_table_json(table, header):
    result = dict(header=header, data=dict())
    for col in table.colnames:
        result['data'][col] = table[col].tolist()

    return result

def render_table_ascii(table, ascii_format):
    with io.StringIO() as buffer:
        table.write(buffer, format=ascii_format)
        result = buffer.getvalue()

    response = make_response(result)

    user_agent = request.headers.get('User-Agent', '').lower()
    ### print(user_agent)
    if 'wget' in user_agent or 'curl' in user_agent:
        if ascii_format == 'csv':
            filename = 'desi-targets.csv'
        else:
            filename = 'desi-targets.txt'

        response.headers['Content-Disposition'] = f'attachment; filename={filename}'
        response.headers['Content-Type'] = 'text/plain'
    else:
        response.headers['Content-Type'] = 'text/plain'

    return response

def render_table_fits(table):
    with io.BytesIO() as buffer:
        table.write(buffer, format='fits')
        buffer.seek(0)
        result = buffer.getvalue()

    response = make_response(result)

    response.headers['Content-Disposition'] = f'attachment; filename=desi-targets.fits'
    response.headers['Content-Type'] = 'application/fits'

    return response

def validate_format(format_type, valid_formats):
    """
    Standardize error message if format_type is not in valid_formats

    Args:
        format_type (str): format to be checked; if None, get from URL query options
        valid_formats: list or tuple of valid format strings

    Returns validated format_type, otherwise raise ValueError
    """
    if format_type is None:
        format_type = request.args.get('format', default='html')

    format_type = format_type.lower()
    if format_type not in valid_formats:
        raise ValueError(f"Unsupported table format='{format_type}'; supported formats are {valid_formats}")

    return format_type

def get_table_format():
    """Return format option from URL, while checking that it is valid for tables"""
    return validate_format(format_type=None, valid_formats=('html', 'json', 'fits', 'csv', 'ascii'))

def get_spectra_format():
    """Return format option from URL, while checking that it is valid for spectra"""
    return validate_format(format_type=None, valid_formats=('html', 'fits'))

def _format_radec(ra,dec):
    rastr = f'{ra:.4f}'.rstrip('0').rstrip('.')
    decstr = f'{dec:.4f}'.rstrip('0').rstrip('.')
    return f'({rastr},{decstr})'

def render_table(table, format_type):

    if format_type == 'html':
        specprod = table.meta['SPECPROD']
        header = f'DESI {specprod} production'
        if ('RA' in table.meta) and ('DEC' in table.meta) and ('RADIUS' in table.meta):
            ra = table.meta['RA']
            dec = table.meta['DEC']
            radius = table.meta['RADIUS']
            radius_str = f'{radius:.1f}'.rstrip('0').rstrip('.')  # drop trailing .0
            description = f'{len(table)} targets within {radius_str} arcsec of RA,dec={_format_radec(ra,dec)}'
        else:
            description = f'{len(table)} targets'

        return render_table_html(table, header, description)

    elif format_type == 'json':
        return render_table_json(table, table.meta)

    elif format_type in ('csv', 'ascii'):
        return render_table_ascii(table, format_type)

    elif format_type in ('fits'):
        return render_table_fits(table)

    else:
        #- Upstream code should have validated format_type, but catch here just in case
        msg = f"Unsupported format='{format_type}'"
        return render_template("error.html", code=400, summary='Bad Request', message=msg), 400

def get_filters():
    """
    Return URL filter args in structure for inspector.io.filter_table
    """
    return request.args.to_dict(flat=False)


#-------------------------------------------------------------------------
#- Inventory of targets

def render_targets(specprod, specgroup, radec=None, targetids=None):
    """
    required: specprod, specgroup; plus radec OR targetids (but not both)
    """
    try:
        format_type = get_table_format()
        filters = get_filters()
        t = load_targets(specprod, specgroup, radec=radec, targetids=targetids, filters=filters)
    except ValueError as err:
        return render_template("error.html", code=400, summary='Bad Request', message=str(err)), 400

    return render_table(t, format_type)

@app.route("/<string:specprod>/targets/radec/<string:radec>")
@app.route("/<string:specprod>/targets/healpix/radec/<string:radec>")
@conditional_auth
def target_healpix_radec(specprod, radec):
    return render_targets(specprod, specgroup='healpix', radec=radec)

@app.route("/<string:specprod>/targets/<string:targetids>")
@app.route("/<string:specprod>/targets/healpix/<string:targetids>")
@conditional_auth
def targets_healpix_targetids(specprod, targetids):
    return render_targets(specprod, specgroup='healpix', targetids=targetids)

@app.route("/<string:specprod>/targets/tiles/radec/<string:radec>")
@conditional_auth
def targets_tiles_radec(specprod, radec):
    return render_targets(specprod, specgroup='tiles', radec=radec)

@app.route("/<string:specprod>/targets/tiles/<string:targetids>")
@conditional_auth
def targets_tiles_targetids(specprod, targetids):
    return render_targets(specprod, specgroup='tiles', targetids=targetids)

@app.route("/<string:specprod>/targets/<int:tileid>/<string:fibers>")
@app.route("/<string:specprod>/targets/tiles/<int:tileid>/<string:fibers>")
@conditional_auth
def targets_tiles_fibers(specprod, tileid, fibers):
    specprod = standardize_specprod(specprod)
    try:
        format_type = get_table_format()
        fibers = parse_fibers(fibers)
    except ValueError as err:
        return render_template("error.html", code=400, summary='Bad Request', message=str(err)), 400

    if np.min(fibers)<0 or np.max(fibers)>=5000:
        msg = 'Fibers must be in 0 <= FIBER < 5000'
        return render_template("error.html", code=400, summary='Bad Request', message=msg), 400

    #- Find LASTNIGHT for this tile
    try:
        lastnight = get_lastnight(tileid, specprod=specprod)
    except ValueError as err:
        msg = f'Tile {tileid} not found in {specprod} production'
        return render_template("error.html", code=404, summary='Not Found', message=str(err)), 404

    #- read fibermaps to find the TARGETIDs for these fibers
    fiber2targetid = dict()
    for petal in np.unique(np.asarray(fibers)//500):
        coaddfile = findfile('coadd', tile=tileid, night=lastnight, groupname='cumulative', spectrograph=petal, specprod=specprod)
        fm = fitsio.read(coaddfile, 'FIBERMAP', columns=('TARGETID', 'FIBER'))
        keep = np.isin(fm['FIBER'], fibers)
        for tid, fiber in zip(fm['TARGETID'][keep], fm['FIBER'][keep]):
            fiber2targetid[fiber] = tid

    #- Create table to render
    targetcat = Table()
    targetcat.meta['SPECPROD'] = specprod
    targetcat['TARGETID'] = [fiber2targetid[f] for f in fibers]
    targetcat['TILEID'] = tileid
    targetcat['LASTNIGHT'] = lastnight
    targetcat['FIBER'] = fibers

    filters = get_filters()
    targetcat = filter_table(add_zcat_columns(targetcat, specprod), filters=filters)

    return render_table(targetcat, format_type)

#-------------------------------------------------------------------------
#- Spectra

def render_spectra_plot(spectra):
    with tempfile.TemporaryDirectory() as tmpdir:

        if 'description' in spectra.meta:
            title = spectra.meta['description']
        else:
            title = 'DESI spectra'

        if 'plotnoise' in request.args and request.args['plotnoise'].lower() in ('', '1', 'true'):
            with_noise = True
        else:
            with_noise = False

        prospectfile = f'{tmpdir}/prospect.html'
        print(f'Writing prospect spectra to {prospectfile}')

        print(f'{len(spectra.fibermap)=}  {len(spectra.redshifts)=}')

        plotspectra(
                spectra,
                zcatalog=spectra.redshifts,
                outfile=prospectfile,
                title=title,
                with_vi_widgets=False,
                with_full_2ndfit=False,
                with_thumb_tab=False,
                with_coaddcam=False,
                with_noise=with_noise,
                num_approx_fits=0,
                top_metadata=['TARGETID', 'TILEID', 'FIBER']
            )

        with open(prospectfile, 'r') as fp:
            html_content = fp.read()

    download_url = _current_url_as_format('fits')
    table_url = _current_url_as_format('html').replace('/spectra/', '/targets/')

    footer = f"""
<span style="font-family: Helvetica, Arial, sans-serif; font-size: 14px;">
<hr/>
DESI Data Inspector:
    <a href='/'>Home</a>
    <a href="{download_url}">Download spectra</a>
    <a href="{table_url}">View target table</a>
</span>
"""

    html_content = html_content.replace('</body>', f'{footer}</body>')

    return html_content, 200, {'Content-Type': 'text/html'}

def render_spectra_fits(spectra):
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = f'{tmpdir}/spectra.fits'
        write_spectra(filename, spectra)

        with open(filename, 'rb') as fp:
            data = fp.read()

    #- make a probably unique filename (1/million chance of collision from 100 downloads)
    blat = hashlib.sha256(request.url.encode()).hexdigest()[0:8]
    filename = f'desi-spectra-{blat}.fits'

    response = make_response(data)
    response.headers['Content-Disposition'] = f'attachment; filename={filename}'
    response.headers['Content-Type'] = 'application/fits'

    return response

def render_spectra(specprod, specgroup, radec=None, targetids=None):
    try:
        format_type = get_spectra_format()
        spectra = load_spectra(specprod, specgroup=specgroup, radec=radec, targetids=targetids)
    except ValueError as err:
        return render_template("error.html", code=400, summary='Bad Request', message=str(err)), 400

    #- Handle case of no spectra found
    if spectra is None:
        if radec is not None:
            ra,dec,radius = validate_radec(radec)
            msg = f'No targets found within {radius} arcsec of RA,dec=({ra},{dec})'
        else:
            msg = f'No targets found with TARGETIDs={targetids}'

        return render_template("error.html", code=400, summary='Not Found', message=str(msg)), 404

    if format_type == 'html':
        if radec is not None:
            ra,dec,radius = validate_radec(radec)
            spectra.meta['description'] = f'{len(spectra)} targets within {radius:.1f} arcsec of RA,dec=({ra:.4f},{dec:.4f})'

        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        #- Upstream code should have validated format_type, but catch here just in case
        msg = f"Unsupportedformat='{format_type}'"
        return render_template("error.html", code=400, summary='Bad Request', message=msg), 400


@app.route("/<string:specprod>/spectra/radec/<string:radec>")
@app.route("/<string:specprod>/spectra/healpix/radec/<string:radec>")
@conditional_auth
def spectra_healpix_radec(specprod, radec):
    return render_spectra(specprod, specgroup='healpix', radec=radec)

@app.route("/<string:specprod>/spectra/<string:targetids>")
@app.route("/<string:specprod>/spectra/healpix/<string:targetids>")
@conditional_auth
def spectra_healpix_targetids(specprod, targetids):
    return render_spectra(specprod, specgroup='healpix', targetids=targetids)

@app.route("/<string:specprod>/spectra/tiles/radec/<string:radec>")
@conditional_auth
def spectra_tiles_radec(specprod, radec):
    return render_spectra(specprod, specgroup='tiles', radec=radec)

@app.route("/<string:specprod>/spectra/tiles/<string:targetids>")
@conditional_auth
def spectra_tiles_targetids(specprod, targetids):
    return render_spectra(specprod, specgroup='tiles', targetids=targetids)

@app.route("/<string:specprod>/spectra/<int:tileid>/<string:fibers>")
@app.route("/<string:specprod>/spectra/tiles/<int:tileid>/<string:fibers>")
@conditional_auth
def spectra_tiles_fibers(specprod, tileid, fibers):
    specprod = standardize_specprod(specprod)
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return render_template("error.html", code=400, summary='Not Found', message=str(err)), 400

    fibers = parse_fibers(fibers)
    if len(fibers) > MAX_SPECTRA:
        msg = MAX_SPECTRA_ERROR_MESSAGE.format(len(fibers), MAX_SPECTRA)
        return render_template("error.html", code=400, summary='Bad Request', message=msg), 400

    #- Find LASTNIGHT for this tile
    try:
        lastnight = get_lastnight(tileid, specprod=specprod)
    except ValueError:
        msg = f'Tile {tileid} not found in {specprod} production'
        return render_template("error.html", code=404, summary='Not Found', message=msg), 404

    targetcat = Table()
    targetcat['FIBER'] = fibers
    targetcat['LASTNIGHT'] = lastnight
    targetcat['TILEID'] = tileid

    filters = get_filters()
    targetcat = filter_table(add_zcat_columns(targetcat, specprod), filters)

    print(f'Reading {len(targetcat)} spectra')
    print(targetcat)
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        #- Upstream code should have validated format_type, but catch here just in case
        msg = f'Unrecognized format {format_type}'
        return render_template("error.html", code=400, summary='Bad Request', message=msg), 400

if __name__ == "__main__":
    app.run(debug=True, port=5500)






