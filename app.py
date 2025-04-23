import io
import os
import glob
import tempfile
from urllib.parse import urlencode
import hashlib
import uuid

from flask import Flask, request, jsonify, render_template, make_response
from astropy.table import Table

from desispec import inventory
from desispec.io import read_spectra_parallel, write_spectra, specprod_root

from prospect.viewer import plotspectra

app = Flask(__name__)
app.url_map.strict_slashes = False

MAX_RADIUS=1800  # 0.5 deg
MAX_RADIUS_ERROR_MESSAGE = f'Please limit your search to radius &le; {MAX_RADIUS} arcsec'

MAX_SPECTRA=1000
MAX_SPECTRA_ERROR_MESSAGE = '{} spectra is more than we can realistically display; please limit your search to fewer than {} spectra'

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
    root_url = request.root_url.rstrip('/')  #- App URL without subpages or query options
    with io.StringIO() as buffer:
        table.write(buffer, format='ascii.html')
        table_html = buffer.getvalue()

    ntargets = len(table)
    nmax = 10000
    if ntargets>nmax:
        msg = f"Your query resulted in {ntargets} targets, which is too many for us to realistically display."
        msg += f" Please limit your query to return less than {nmax} targets."
        return msg

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
    print(user_agent)
    if 'wget' in user_agent or 'curl' in user_agent:
        response.headers['Content-Disposition'] = f'attachment; filename=desi-inventory.txt'
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

    response.headers['Content-Disposition'] = f'attachment; filename=desi-inventory.fits'
    response.headers['Content-Type'] = 'application/fits'

    return response

def get_table_format():
    format_type = request.args.get('format', default='html')
    valid_formats = ('html', 'json', 'fits', 'csv', 'ascii')
    if format_type not in valid_formats:
        raise ValueError(f"Unsupported format='{format_type}'; supported formats are {valid_formats}")

    return format_type

def get_spectra_format():
    format_type = request.args.get('format', default='html')
    valid_formats = ('html', 'fits')
    if format_type not in valid_formats:
        raise ValueError(f"Unsupported format='{format_type}'; supported formats are {valid_formats}")

    return format_type

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
        return f"Unsupported format='{format_type}'", 400

#-------------------------------------------------------------------------
#- Inventory of targets

@app.route("/<string:specprod>/targets/radec/<string:radec>")
@app.route("/<string:specprod>/targets/healpix/radec/<string:radec>")
def healpix_radec_targets(specprod, radec):
    try:
        format_type = get_table_format()
    except ValueError as err:
        return str(err), 400

    ra,dec,radius = inventory.parse_radec_string(radec)
    if radius > MAX_RADIUS:
        return MAX_RADIUS_ERROR_MESSAGE

    t = inventory.target_healpix(radec=(ra,dec,radius), specprod=specprod)
    return render_table(t, format_type)


@app.route("/<string:specprod>/targets/<string:targetids>")
@app.route("/<string:specprod>/targets/healpix/<string:targetids>")
def healpix_targets(specprod, targetids):
    try:
        format_type = get_table_format()
    except ValueError as err:
        return str(err), 400

    targetids = list(map(int, targetids.split(',')))
    t = inventory.target_healpix(targetids=targetids, specprod=specprod)
    return render_table(t, format_type)


@app.route("/<string:specprod>/targets/tiles/radec/<string:radec>")
def tiles_radec_targets(specprod, radec):
    try:
        format_type = get_table_format()
    except ValueError as err:
        return str(err), 400

    ra,dec,radius = inventory.parse_radec_string(radec)
    if radius > MAX_RADIUS:
        return MAX_RADIUS_ERROR_MESSAGE

    t = inventory.target_tiles(radec=(ra,dec,radius), specprod=specprod)
    return render_table(t, format_type)


@app.route("/<string:specprod>/targets/tiles/<string:targetids>")
def tiles_targets(specprod, targetids):
    try:
        format_type = get_table_format()
    except ValueError as err:
        return str(err), 400

    targetids = list(map(int, targetids.split(',')))
    t = inventory.target_tiles(targetids=targetids, specprod=specprod)

    return render_table(t, format_type)

#-------------------------------------------------------------------------
#- Spectra

def render_spectra_plot(spectra):
    with tempfile.TemporaryDirectory() as tmpdir:

        if 'description' in spectra.meta:
            title = spectra.meta['description']
        else:
            title = 'DESI spectra'

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
                num_approx_fits=0,
            )

        with open(prospectfile, 'r') as fp:
            html_content = fp.read()

    download_url = _current_url_as_format('fits')
    table_url = _current_url_as_format('html').replace('/spectra/', '/targets/')

    html_content = html_content.replace('</body>',
                                        f'<hr/><a href="{download_url}">Download spectra</a>; <a href="{table_url}">View target table</a>\n</body>')

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

@app.route("/<string:specprod>/spectra/radec/<string:radec>")
@app.route("/<string:specprod>/spectra/healpix/radec/<string:radec>")
def plot_healpix_spectra_radec(specprod, radec):
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return str(err), 400

    ra,dec,radius = inventory.parse_radec_string(radec)
    targetcat = inventory.target_healpix(radec=(ra,dec,radius), specprod=specprod)

    num_spectra = len(targetcat)
    if num_spectra == 0:
        return f'No targets found within {radius} arcsec of RA,dec=({ra},{dec})', 404
    elif num_spectra > MAX_SPECTRA:
        return MAX_SPECTRA_ERROR_MESSAGE.format(num_spectra, MAX_SPECTRA)

    print(f'Reading {num_spectra} spectra')
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        spectra.meta['description'] = f'{len(spectra)} targets within {radius:.1f} arcsec of RA,dec=({ra:.4f},{dec:.4f})'
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        return f'Unrecognized format {format_type}', 400


@app.route("/<string:specprod>/spectra/<string:targetids>")
@app.route("/<string:specprod>/spectra/healpix/<string:targetids>")
def plot_healpix_spectra_targetids(specprod, targetids):
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return str(err), 400

    targetids = list(map(int, targetids.split(',')))
    targetcat = inventory.target_healpix(targetids=targetids, specprod=specprod)

    num_spectra = len(targetcat)
    if num_spectra == 0:
        return f'No targets found within {radius} arcsec of RA,dec=({ra},{dec})', 404
    elif num_spectra > MAX_SPECTRA:
        return MAX_SPECTRA_ERROR_MESSAGE.format(num_spectra, MAX_SPECTRA)

    print(f'Reading {num_spectra} spectra')
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        return f'Unrecognized format {format_type}', 400


@app.route("/<string:specprod>/spectra/tiles/radec/<string:radec>")
def plot_tiles_spectra_radec(specprod, radec):
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return str(err), 400

    ra,dec,radius = inventory.parse_radec_string(radec)
    targetcat = inventory.target_tiles(radec=(ra,dec,radius), specprod=specprod)

    num_spectra = len(targetcat)
    if num_spectra == 0:
        return f'No targets found within {radius} arcsec of RA,dec=({ra},{dec})', 404
    elif num_spectra > MAX_SPECTRA:
        return MAX_SPECTRA_ERROR_MESSAGE.format(num_spectra, MAX_SPECTRA)

    print(f'Reading {num_spectra} spectra')
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        spectra.meta['description'] = f'{len(spectra)} targets within {radius:.1f} arcsec of RA,dec=({ra:.4f},{dec:.4f})'
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        return f'Unrecognized format {format_type}', 400


@app.route("/<string:specprod>/spectra/tiles/<string:targetids>")
def plot_tiles_spectra_targetids(specprod, targetids):
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return str(err), 400

    targetids = list(map(int, targetids.split(',')))
    targetcat = inventory.target_tiles(targetids=targetids, specprod=specprod)

    num_spectra = len(targetcat)
    if num_spectra == 0:
        return f'No targets found with TARGETIDs={targetids}'
    elif num_spectra > MAX_SPECTRA:
        return MAX_SPECTRA_ERROR_MESSAGE.format(num_spectra, MAX_SPECTRA)

    print(f'Reading {num_spectra} spectra')
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        return f'Unrecognized format {format_type}', 400

@app.route("/<string:specprod>/spectra/tiles/<int:tileid>/<string:fibers>")
def plot_tiles_spectra_fibers(specprod, tileid, fibers):
    try:
        format_type = get_spectra_format()
    except ValueError as err:
        return str(err), 400

    fibers = list(map(int, fibers.split(',')))
    if len(fibers) > MAX_SPECTRA:
        return MAX_SPECTRA_ERROR_MESSAGE.format(len(fibers), MAX_SPECTRA)

    #- presumably we won't be running this code past the year 2100
    nightdirs = sorted(glob.glob(specprod_root(specprod)+f'/tiles/cumulative/{tileid}/20??????'))
    lastnight = int(os.path.basename(nightdirs[-1]))

    targetcat = Table()
    targetcat['FIBER'] = fibers
    targetcat['LASTNIGHT'] = lastnight
    targetcat['TILEID'] = tileid

    print(f'Reading {len(targetcat)} spectra')
    print(targetcat)
    spectra = read_spectra_parallel(targetcat, specprod=specprod, rdspec_kwargs=dict(return_redshifts=True))

    if format_type == 'html':
        return render_spectra_plot(spectra)
    elif format_type == 'fits':
        return render_spectra_fits(spectra)
    else:
        return f'Unrecognized format {format_type}', 400

if __name__ == "__main__":
    app.run(debug=True)






