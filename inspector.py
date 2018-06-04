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
from bokeh.models import Legend, Range1d
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.layouts import row, column
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Div
import bokeh.palettes
import bokeh.events

import desispec.io
from desispec.spectra import Spectra
from desispec.interpolation import resample_flux
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask

import redrock.templates
from redrock.rebin import trapz_rebin

def read_templates():
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


def load_spectra(specfile, zbestfile=None):
    '''TODO: document'''
    if zbestfile is None:
        specdir, basename = os.path.split(specfile)
        if not basename.startswith('spectra'):
            raise ValueError("Can't derive zbest filename if spectra filename {} doesn't match spectra*.fits".format(basename))
        zbestfile = os.path.join(specdir, basename.replace('spectra', 'zbest'))

        zbest = Table.read(zbestfile, 'ZBEST')
        spectra = coadd_targets(desispec.io.read_spectra(specfile),
                targetids=zbest['TARGETID'])

    return Inspector(spectra, zbest)

class Inspector(): 
    def __init__(self, spectra, zbest):
        '''TODO: document'''
        self.zbest = zbest
        self.spectra = spectra
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
        self.print_targets_info()
        output_notebook()

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
            
        p.yaxis.axis_label = 'Flux [1e-17 erg/s/cm2/A]'
        p.xaxis.axis_label = 'Wavelength [A]'
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
                y_range=p.y_range, output_backend="webgl", tools=[])
        for channel in ['b', 'r', 'z']:
            pz.line('wave', 'flux', source=self.data[channel],
                    line_color=colors[channel], line_width=1, line_alpha=1.0)
            pz.line('wave', 'model', source=self.data[channel],
                    line_color='black', line_width=1, alpha=1.0)

        z = self.zbest['Z'][self.izbest]
        pz.x_range.start = 3727*(1+z) - 100
        pz.x_range.end = 3727*(1+z) + 100
        self.pz = pz

        #- Callback to update zoom window x-range 
        def callback(zoomplot):
            return CustomJS(args=dict(xr=zoomplot.x_range), code="""
                xr.start = cb_obj.x - 100;
                xr.end = cb_obj.x + 100;
            """)
    
        p.js_on_event(bokeh.events.MouseMove, callback(pz))

        #- Text area
        self.infodiv = Div(text='Hello<br/>There')

        self.plot_handle = show(row(p, column(pz, widgetbox(self.infodiv))),
                notebook_handle=True)
        # self.plot_handle = show(p, notebook_handle=True)

        #- Repeatd code...
        self._update()

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
        
        z = zb['Z']
        self.pz.x_range.start = 3727*(1+z) - 100
        self.pz.x_range.end = 3727*(1+z) + 100

        fibermap = self.spectra.fibermap[self.ispec]

        title = '{} z={:.4f} zwarn={}'.format(
            zb['SPECTYPE'], zb['Z'], zb['ZWARN'])
        self.p.title.text = title

        info = list()
        info.append('{}/{}'.format(self.ispec+1, self.spectra.num_spectra()))
        info.append('ID {}'.format(zb['TARGETID']))
        if zb['ZWARN']:
            info.append('<font color="orangered"><b>ZWARN={}</b></font>'.format(
                zb['ZWARN']))
        info.append('<b>DESI_TARGET</b><br/>&nbsp;&nbsp;{}'.format(
            ' '.join(desi_mask.names(fibermap['DESI_TARGET']))))
        if fibermap['BGS_TARGET']:
            info.append('<b>BGS_TARGET</b><br/>&nbsp;&nbsp;{}'.format(
                ' '.join(bgs_mask.names(fibermap['BGS_TARGET']))))

        if fibermap['MWS_TARGET']:
            info.append('<b>MWS_TARGET</b><br/>&nbsp;&nbsp;{}'.format(
                ' '.join(mws_mask.names(fibermap['MWS_TARGET']))))

        info = '<br/>\n'.join(info)
        self.infodiv.text = '''<p style="font-size:0.8em">{}</p>'''.format(info)

        push_notebook(handle=self.plot_handle)


