'''
Inspector: it's for inspecting DESI spectra

Stephen Bailey
Spring 2018
'''

import os, sys, glob
import numpy as np

from astropy.table import Table

from bokeh.io import push_notebook, show, output_notebook
from bokeh.plotting import figure
from bokeh.models import (CustomJS, ColumnDataSource, Label, Legend,
                          Range1d, Slider, Span, Arrow, VeeHead)
from bokeh.layouts import row, column, widgetbox
from bokeh.models.widgets import Div
import bokeh.palettes
import bokeh.events

import desispec.io
from desispec.spectra import Spectra
from desispec.interpolation import resample_flux
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

import redrock.templates
from redrock.rebin import trapz_rebin


lines = [
    #
    # This is the set of emission lines from the spZline files.
    # See $IDLSPEC2D_DIR/etc/emlines.par
    # Wavelengths are in air for lambda > 2000, vacuum for lambda < 2000.
    #
    {"name" : "Lyα",      "longname" : "Lyman α",        "lambda" : 1215.67,  "emission": True },
    {"name" : "N V",      "longname" : "N V 1240",       "lambda" : 1240.81,  "emission": True },
    {"name" : "C IV",     "longname" : "C IV 1549",      "lambda" : 1549.48,  "emission": True },
    {"name" : "He II",    "longname" : "He II 1640",     "lambda" : 1640.42,  "emission": True },
    {"name" : "C III]",   "longname" : "C III] 1908",    "lambda" : 1908.734, "emission": True },
    {"name" : "Mg II",    "longname" : "Mg II 2799",     "lambda" : 2799.49,  "emission": True },
    {"name" : "[O II]",   "longname" : "[O II] 3725",    "lambda" : 3726.032, "emission": True },
    {"name" : "[O II]",   "longname" : "[O II] 3727",    "lambda" : 3728.815, "emission": True },
    {"name" : "[Ne III]", "longname" : "[Ne III] 3868",  "lambda" : 3868.76,  "emission": True },
    {"name" : "Hζ",       "longname" : "Balmer ζ",       "lambda" : 3889.049, "emission": True },
    {"name" : "[Ne III]", "longname" : "[Ne III] 3970",  "lambda" : 3970.00,  "emission": True },
    {"name" : "Hε",       "longname" : "Balmer ε",       "lambda" : 3970.072, "emission": True },
    {"name" : "Hδ",       "longname" : "Balmer δ",       "lambda" : 4101.734, "emission": True },
    {"name" : "Hγ",       "longname" : "Balmer γ",       "lambda" : 4340.464, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 4363",   "lambda" : 4363.209, "emission": True },
    {"name" : "He II",    "longname" : "He II 4685",     "lambda" : 4685.68,  "emission": True },
    {"name" : "Hβ",       "longname" : "Balmer β",       "lambda" : 4861.325, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 4959",   "lambda" : 4958.911, "emission": True },
    {"name" : "[O III]",  "longname" : "[O III] 5007",   "lambda" : 5006.843, "emission": True },
    {"name" : "He II",    "longname" : "He II 5411",     "lambda" : 5411.52,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 5577",     "lambda" : 5577.339, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 5755",    "lambda" : 5754.59,  "emission": True },
    {"name" : "He I",     "longname" : "He I 5876",      "lambda" : 5875.68,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 6300",     "lambda" : 6300.304, "emission": True },
    {"name" : "[S III]",  "longname" : "[S III] 6312",   "lambda" : 6312.06,  "emission": True },
    {"name" : "[O I]",    "longname" : "[O I] 6363",     "lambda" : 6363.776, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 6548",    "lambda" : 6548.05,  "emission": True },
    {"name" : "Hα",       "longname" : "Balmer α",       "lambda" : 6562.801, "emission": True },
    {"name" : "[N II]",   "longname" : "[N II] 6583",    "lambda" : 6583.45,  "emission": True },
    {"name" : "[S II]",   "longname" : "[S II] 6716",    "lambda" : 6716.44,  "emission": True },
    {"name" : "[S II]",   "longname" : "[S II] 6730",    "lambda" : 6730.82,  "emission": True },
    {"name" : "[Ar III]", "longname" : "[Ar III] 7135",  "lambda" : 7135.790, "emission": True },
    #
    # Absorption lines
    #
    {"name" : "Hζ",   "longname" : "Balmer ζ",         "lambda" : 3889.049, "emission": False },
    {"name" : "K",    "longname" : "K (Ca II 3933)",   "lambda" : 3933.7,   "emission": False },
    {"name" : "H",    "longname" : "H (Ca II 3968)",   "lambda" : 3968.5,   "emission": False },
    {"name" : "Hε",   "longname" : "Balmer ε",         "lambda" : 3970.072, "emission": False },
    {"name" : "Hδ",   "longname" : "Balmer δ",         "lambda" : 4101.734, "emission": False },
    {"name" : "G",    "longname" : "G (Ca I 4307)",    "lambda" : 4307.74,  "emission": False },
    {"name" : "Hγ",   "longname" : "Balmer γ",         "lambda" : 4340.464, "emission": False },
    {"name" : "Hβ",   "longname" : "Balmer β",         "lambda" : 4861.325, "emission": False },
    {"name" : "Mg I", "longname" : "Mg I 5175",        "lambda" : 5175.0,   "emission": False },
    {"name" : "D2",   "longname" : "D2 (Na I 5889)",   "lambda" : 5889.95,  "emission": False },
    # {"name" : "D",    "longname" : "D (Na I doublet)", "lambda": 5892.9,   "emission": False },
    {"name" : "D1",   "longname" : "D1 (Na I 5895)",   "lambda" : 5895.92,  "emission": False },
    {"name" : "Hα",   "longname" : "Balmer α",         "lambda" : 6562.801, "emission": False },
    ]

def airtovac(l):
    """Convert air wavelengths to vacuum wavelengths. Don't convert less than 2000 Å.

    Parameters
    ----------
    l : :class:`float`
        Wavelength [Å] of the line in air.

    Returns
    -------
    :class:`float`
        Wavelength [Å] of the line in vacuum.
    """
    if l < 2000.0:
        return l;
    vac = l
    for iter in range(2):
        sigma2 = (1.0e4/vac)*(1.0e4/vac)
        fact = 1.0 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
        vac = l*fact
    return vac

def read_templates():
    """Retrieve redrock templates.

    Returns
    -------
    :class:`dict`
        A dictionary keyed on (type, subtype).
    """
    #- redirect stdout to silence chatty redrock
    saved_stdout = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    try:
        templates = dict()
        for filename in redrock.templates.find_templates():
            t = redrock.templates.Template(filename)
            templates[(t.template_type, t.sub_type)] = t
    except Exception as err:
        sys.stdout = saved_stdout
        raise(err)

    sys.stdout = saved_stdout
    return templates

def coadd(wave, flux, ivar, rdat):
    nspec, nwave = flux.shape
    unweightedflux = np.zeros(nwave, dtype=flux.dtype)
    weightedflux = np.zeros(nwave, dtype=flux.dtype)
    weights = np.zeros(nwave, dtype=flux.dtype)
    outrdat = np.zeros(rdat[0].shape, dtype=rdat.dtype)
    for i in range(nspec):
        unweightedflux += flux[i]
        weightedflux += flux[i] * ivar[i]
        weights += ivar[i]
        outrdat += rdat[i] * ivar[i]

    isbad = (weights == 0)
    outflux = weightedflux / (weights + isbad)
    outflux[isbad] = unweightedflux[isbad] / nspec

    outrdat /= (weights + isbad)
    outivar = weights

    return wave, outflux, outivar, outrdat

def coadd_targets(spectra, targetids=None):
    '''
    Coadds individual exposures of the same targets; returns new Spectra object
    '''
    if targetids is None:
        targetids = spectra.target_ids()

    ntargets = spectra.num_targets()
    wave = dict()
    flux = dict()
    ivar = dict()
    rdat = dict()
    for channel in spectra.bands:
        wave[channel] = spectra.wave[channel].copy()
        nwave = len(wave[channel])
        flux[channel] = np.zeros((ntargets, nwave))
        ivar[channel] = np.zeros((ntargets, nwave))
        ndiag = spectra.resolution_data[channel].shape[1]
        rdat[channel] = np.zeros((ntargets, ndiag, nwave))

    fibermap = Table(dtype=spectra.fibermap.dtype)
    for i, targetid in enumerate(targetids):
        ii = np.where(spectra.fibermap['TARGETID'] == targetid)[0]
        fibermap.add_row(spectra.fibermap[ii[0]])
        for channel in spectra.bands:
            if len(ii) > 1:
                outwave, outflux, outivar, outrdat = coadd(
                    spectra.wave[channel],
                    spectra.flux[channel][ii],
                    spectra.ivar[channel][ii],
                    spectra.resolution_data[channel][ii]
                    )
            else:
                outwave, outflux, outivar, outrdat = (
                    spectra.wave[channel],
                    spectra.flux[channel][ii[0]],
                    spectra.ivar[channel][ii[0]],
                    spectra.resolution_data[channel][ii[0]]
                    )

            flux[channel][i] = outflux
            ivar[channel][i] = outivar
            rdat[channel][i] = outrdat

    return Spectra(spectra.bands, wave, flux, ivar,
            mask=None, resolution_data=rdat, fibermap=fibermap)


class Inspector(object):
    """An interface to plotting spectra with Bokeh.

    Parameters
    ----------
    basedir : :class:`str`
        Path to the directory containing spectrum files.
    """

    def __init__(self, basedir):
        specfile = glob.glob(os.path.join(basedir, 'spectra-*.fits'))[0]
        zbestfile = glob.glob(os.path.join(basedir, 'zbest-*.fits'))[0]
        self.zbest = Table.read(zbestfile, 'ZBEST')
        self.spectra = coadd_targets(desispec.io.read_spectra(specfile),
                targetids=self.zbest['TARGETID'])
        self.templates = read_templates()
        self.nspec = len(self.zbest)

        sp = self.spectra
        assert np.all(sp.target_ids() == self.zbest['TARGETID'])
        assert np.all(sp.target_ids() == sp.fibermap['TARGETID'])

        self.data = dict()     #- high resolution
        self.xdata = dict()    #- low resolution
        for channel in ['b', 'r', 'z']:
            izbest, (xwave, xflux, xmodel), (wave, flux, model) = \
                    self._get_data(0, channel)
            self.data[channel] = ColumnDataSource(
                    dict(wave=wave, flux=flux, model=model))
            self.xdata[channel] = ColumnDataSource(
                    dict(wave=xwave, flux=xflux, model=xmodel))

        self.ispec = 0
        self.izbest = izbest
        self._emission = False
        self._absorption = False
        self.print_targets_info()
        output_notebook()

    @property
    def z(self):
        """This is used in several places, so best to make it a property.
        """
        return self.zbest[self.izbest]['Z']

    def select(self, targetids):
        ii = np.in1d(self.zbest['TARGETID'], targetids)
        self.zbest = self.zbest[ii]
        self.spectra = self.spectra.select(targets=targetids)
        self.nspec = len(self.zbest)
        self.ispec = 0
        self.izbest = 0
        self.print_targets_info()

    def print_targets_info(self):
        ntargets = self.spectra.num_targets()
        fm = self.spectra.fibermap
        nexp = len(fm)
        nelg = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.ELG)
        nlrg = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.LRG)
        nqso = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.QSO)
        nbgs = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.BGS_ANY)
        nmws = np.count_nonzero(fm['DESI_TARGET'] & desi_mask.MWS_ANY)
        print('{} targets'.format(ntargets), end='')
        print(' including {} ELG, {} LRG, {} QSO, {} BGS, {} MWS'.format(
            nelg, nlrg, nqso, nbgs, nmws))

    def plot(self):
        tools = 'pan,box_zoom,wheel_zoom,undo,redo,reset,save'
        self.p = p = figure(plot_height=400, plot_width=800,
                        output_backend="webgl",
                        toolbar_location='above', tools=tools)

        p.toolbar.active_drag = p.tools[0]    #- pan
        p.toolbar.active_scroll = p.tools[2]  #- wheel zoom

        colors = dict(b='#1f77b4', r='#d62728', z='maroon')
        flux_lines = list()
        model_lines = list()
        for channel in ['b', 'r', 'z']:
            flux_lines.append(p.line('wave', 'flux',
                source=self.xdata[channel],
                line_color=colors[channel], line_width=1, alpha=1.0))
            model_lines.append(p.line('wave', 'model',
                source=self.xdata[channel],
                line_color='black', line_width=1, alpha=1.0))

        xtmp = [self.xdata['b'].data['wave'][0],
                self.xdata['z'].data['wave'][-1]]
        ytmp = [0,0]
        p.line(xtmp, ytmp, color='black', alpha=0.5)

        p.yaxis.axis_label = 'Flux [10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹]'
        p.xaxis.axis_label = 'Wavelength [Å]'
        p.xaxis.axis_label_text_font_style = 'normal'
        p.yaxis.axis_label_text_font_style = 'normal'
        p.min_border_left = 60
        p.min_border_bottom = 40
        self._set_ylim()

        legend = Legend(items=[
            ("flux",  flux_lines),
            ("model", model_lines),
        ])
        p.add_layout(legend, 'center')
        p.legend.click_policy = 'hide'    #- or 'mute'

        #- Zoom plot
        pz = figure(title=None, plot_height=200, plot_width=200,
                y_range=p.y_range, output_backend="webgl",
                toolbar_location='above', tools=[])
        for channel in ['b', 'r', 'z']:
            pz.line('wave', 'flux', source=self.data[channel],
                    line_color=colors[channel], line_width=1, line_alpha=1.0)
            pz.line('wave', 'model', source=self.data[channel],
                    line_color='black', line_width=1, alpha=1.0)

        pz.x_range.start = 3727*(1+self.z) - 100
        pz.x_range.end = 3727*(1+self.z) + 100
        self.pz = pz

        #- Callback to update zoom window x-range
        def callback(zoomplot):
            return CustomJS(args=dict(xr=zoomplot.x_range), code="""
                xr.start = cb_obj.x - 100;
                xr.end = cb_obj.x + 100;
            """)

        p.js_on_event(bokeh.events.MouseMove, callback(pz))
        #
        # Lines.
        #
        self._update_lines()
        #
        # Targeting information
        #
        self.im = figure(title=None, plot_width=256, plot_height=256,
                         x_range=(0, 256), y_range=(0, 256),
                         output_backend="webgl",
                         toolbar_location='above', tools=[])
        u = self._cutout(self.spectra.fibermap[self.ispec]['RA_TARGET'],
                         self.spectra.fibermap[self.ispec]['DEC_TARGET'])
        self.targetdiv = Div(text='Hello<br/>There')

        #- Text area
        self.infodiv = Div(text='Hello<br/>There')

        self.plot_handle = show(row(column(self.im, widgetbox(self.targetdiv)),
                                    p,
                                    column(pz, widgetbox(self.infodiv))),
                                notebook_handle=True)
        # self.plot_handle = show(p, notebook_handle=True)

        #- Repeatd code...
        self._update()

    def _cutout(self, ra, dec, zoom=12, layer='sdss2'):
        """Image plot centered on `ra`, `dec`.
        """
        u = "http://legacysurvey.org/viewer/jpeg-cutout?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}".format(ra, dec, zoom, layer)
        v = "http://legacysurvey.org/viewer/?ra={0:f}&dec={1:f}&zoom={2:d}&layer={3}".format(ra, dec, zoom, layer)
        self.im.image_url([u], 1, 1, 256, 256, anchor='bottom_left')
        return v

    def _set_ylim(self):
        ymin = ymax = 0.0
        for channel in ['b', 'r', 'z']:
            model = self.data[channel].data['model']
            flux = self.data[channel].data['flux']
            ymax = max(ymax, np.max(model)*1.05)
            ymax = max(ymax, np.percentile(flux, 98))
            ymin = min(ymin, np.percentile(flux, 10))
            ymin = min(0, ymin)

        self.p.y_range.start = ymin
        self.p.y_range.end = ymax

    def _get_data(self, ispec, channel):
        wave = self.spectra.wave[channel]
        flux = self.spectra.flux[channel][ispec]
        targetid = self.spectra.fibermap['TARGETID'][ispec]
        izbest = np.where(self.zbest['TARGETID']==targetid)[0][0]
        zb = self.zbest[izbest]
        z = zb['Z']
        spectype = (zb['SPECTYPE'], zb['SUBTYPE'])
        tx = self.templates[spectype]
        coeff = zb['COEFF'][0:tx.nbasis]
        model = tx.flux.T.dot(coeff).T

        ww = np.arange(wave[0], wave[-1], 3)
        xflux = resample_flux(ww, wave, flux)
        xmodel = resample_flux(ww, tx.wave*(1+z), model)

        model = resample_flux(wave, tx.wave*(1+z), model)

        return izbest, (ww, xflux, xmodel), (wave, flux, model)

    def next(self):
        if self.ispec+1 < self.nspec:
            self.ispec += 1
        else:
            print('end of targets')
        self._update()

    def prev(self):
        if self.ispec > 0:
            self.ispec -= 1
        else:
            print('Already at first target')
        self._update()

    def emission(self, toggle=None):
        """Toggle the display of known emission lines.

        Parameters
        ----------
        toggle : :class:`bool`, optional
            ``True`` and ``False`` turn on and off emission lines,
            respectively.  If not set, the state will be set to the
            opposite of the current state.
        """
        if toggle is None:
            self._emission = not self._emission
        else:
            self._emission = bool(toggle)
        self._update_lines()
        push_notebook(handle=self.plot_handle)

    def absorption(self, toggle=None):
        """Toggle the display of known absorption lines.

        Parameters
        ----------
        toggle : :class:`bool`, optional
            ``True`` and ``False`` turn on and off emission lines,
            respectively.  If not set, the state will be set to the
            opposite of the current state.
        """
        if toggle is None:
            self._absorption = not self._absorption
        else:
            self._absorption = bool(toggle)
        self._update_lines()
        push_notebook(handle=self.plot_handle)

    def _update_lines(self):
        for i, l in enumerate(lines):
            shiftedWave = airtovac(l['lambda'])*(1.0 + self.z)
            visible = (self._line_in_range(shiftedWave) and
                       ((l['emission'] and self._emission) or
                        (self._absorption and not l['emission'])))
            for channel in ('b', 'r', 'z'):
                if self.xdata[channel].data['wave'].min() < l < self.xdata[channel].data['wave'].max():
                    shiftedWave_y = np.interp(shiftedWave,
                                              self.xdata[channel].data['wave'],
                                             self.xdata[channel].data['model'])
                    break
            if l['emission']:
                lc = 'blue'
                y_start = shiftedWave_y*1.2
                y_end = shiftedWave_y*1.2 - 0.5
                # y_start = 150
                # y_end = 140
            else:
                lc = 'red'
                y_start = shiftedWave_y*0.8
                y_end = shiftedWave_y*0.8 + 0.5
                # y_start = 50
                # y_end = 60
            if 'span' in l:
                l['span'].x_start = shiftedWave
                l['span'].y_start = y_start
                l['span'].x_end = shiftedWave
                l['span'].y_start = y_end
                # l['span'].location = shiftedWave
                l['span'].visible = visible
                l['label'].x = shiftedWave
                l['label'].visible = visible
            else:
                l['span'] = Arrow(end=VeeHead(line_color=lc, line_alpha=0.3,
                                              fill_color=lc, fill_alpha=0.3),
                                  line_color=lc, line_width=3, line_alpha=0.3,
                                  x_start=shiftedWave, y_start=y_start,
                                  x_end=shiftedWave, y_end=y_end,
                                  visible=visible)
                # l['span'] = Span(location=shiftedWave, dimension='height',
                #                  line_color=lc, line_dash='solid',
                #                  line_width=3, line_alpha=0.3,
                #                  visible=visible)
                l['label'] = Label(x=shiftedWave, y=y_start,
                                   text=l['name'], text_color=lc, text_alpha=0.3,
                                   visible=visible)
                self.p.add_layout(l['span'])
                self.p.add_layout(l['label'])

    def _line_in_range(self, l):
        """True if a spectral line is within the range of the plot.

        Parameters
        ----------
        l : :class:`float`
            Wavelength [Å] of the line to be tested.

        Returns
        -------
        :class:`bool`
            ``True`` if the line should be plotted.
        """
        return self.xdata['b'].data['wave'].min() < l < self.xdata['z'].data['wave'].max()

    def _update(self):
        for channel in ['b', 'r', 'z']:
            izbest, (xwave, xflux, xmodel), (wave, flux, model) = \
                    self._get_data(self.ispec, channel)
            self.xdata[channel].data['wave'] = xwave
            self.xdata[channel].data['flux'] = xflux
            self.xdata[channel].data['model'] = xmodel
            self.data[channel].data['wave'] = wave
            self.data[channel].data['flux'] = flux
            self.data[channel].data['model'] = model

        self.izbest = izbest
        zb = self.zbest[self.izbest]

        self._set_ylim()

        self.pz.x_range.start = 3727*(1 + self.z) - 100
        self.pz.x_range.end = 3727*(1 + self.z) + 100

        self._update_lines()
        fibermap = self.spectra.fibermap[self.ispec]

        title = '{0} z={1:.4f} zwarn={2}'.format(
            zb['SPECTYPE'], zb['Z'], zb['ZWARN'])
        self.p.title.text = title

        info = ['<dl style="font-size:small;">']
        info.append('<dt>Spectrum</dt><dd>{0:d}/{1:d}</dd>'.format(self.ispec+1, self.spectra.num_spectra()))
        info.append('<dt>ID</dt><dd>{0:d}</dd>'.format(zb['TARGETID']))
        if zb['ZWARN']:
            info.append('<dt style="color:orangered;">ZWARN</dt><dd style="color:orangered;">{0:d}</dd>'.format(
                        zb['ZWARN']))
        info.append('<dt>DESI_TARGET</dt><dd>{0}<dd>'.format(
            ' '.join(desi_mask.names(fibermap['DESI_TARGET']))))
        if fibermap['BGS_TARGET']:
            info.append('<dt>BGS_TARGET</dt><dd>{0}</dd>'.format(
                ' '.join(bgs_mask.names(fibermap['BGS_TARGET']))))

        if fibermap['MWS_TARGET']:
            info.append('<dt>MWS_TARGET</dt><dd>{0}</dd>'.format(
                ' '.join(mws_mask.names(fibermap['MWS_TARGET']))))
        info.append('</dl>')
        self.infodiv.text = ''.join(info)

        targeturl = self._cutout(fibermap['RA_TARGET'], fibermap['DEC_TARGET'])
        targetopen = "window.open('{0}', '_blank');".format(targeturl)
        targetinfo = ['<dl style="font-size:small;">']
        targetinfo.append('<dt>Viewer Link</dt><dd><a onclick="{0}">DECaLS Image Viewer</a></dd>'.format(targetopen))
        targetinfo.append('<dt>RA</dt><dd>{0:.6f}</dd>'.format(fibermap['RA_TARGET']))
        targetinfo.append('<dt>DEC</dt><dd>{0:.6f}</dd>'.format(fibermap['DEC_TARGET']))
        targetinfo.append('</dl>')
        self.targetdiv.text = ''.join(targetinfo)


        push_notebook(handle=self.plot_handle)
