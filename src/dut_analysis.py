#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a single DUT
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
import src.bins as bins
from mod.dut_cuts import DUTCut
from plotting.fit import *
from src.analysis import *
from src.calibration import Calibration
from src.cern_run import CERNRun
from src.converter import Converter
from src.currents import Currents
from src.desy_converter import DESYConverter
from src.desy_run import DESYRun
from src.dummy import Dummy
from src.dut import Plane
from utility.utils import *
from utility.affine_transformations import transform


class DUTAnalysis(Analysis):

    def __init__(self, run_number, dut_number, test_campaign, single_mode=True, verbose=True, test=False):

        Analysis.__init__(self, test_campaign, meta_sub_dir='DUT', verbose=verbose)
        self.print_start(run_number)

        # MAIN
        self.Run = self.run(run_number, dut_number, self.BeamTest.Path, self.Config, single_mode)
        self.DUT = self.Run.DUT
        self.Planes = self.init_planes()
        self.Plane = self.DUT.Plane

        # DATA
        self.Converter = self.converter(self.BeamTest.Path, self.Run.Number, self.Config)
        if test:
            return

        self.Dummy = Dummy(self.BeamTest.Path, self.Converter.NTelPlanes, self.Converter.NDUTPlanes, self.Config)
        self.F = self.load_file(test)
        self.T = False

        # INFO
        self.N = self.n
        self.NEvents = self.F['Tracks']['NTracks'].size
        self.NTracks = self.F['Tracks']['Size'].size
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()
        self.Duration = (self.EndTime - self.StartTime).seconds
        self.Surface = False

        # SUBCLASSES
        self.Calibration = Calibration(self.Run)
        self.Cut = DUTCut(self)
        self.Residuals = self.init_residuals()

        self.Tel = self.init_tel()
        self.REF = self.init_ref()
        self.Tracks = self.init_tracks()

        self.Currents = Currents(self)
        self.Cut.make_additional()
        self.Efficiency = self.init_eff()

    def __repr__(self):
        return f'{super().__repr__()} ({ev2str(self.NEvents)} ev)'

    # ----------------------------------------
    # region INIT
    @property
    def run(self):
        return DESYRun if self.BeamTest.Location == 'DESY' else CERNRun

    @property
    def converter(self):
        return DESYConverter if self.BeamTest.Location == 'DESY' else Converter

    def init_planes(self):
        n_tel, n_dut = [self.Config.get_value(section, 'planes', dtype=int) for section in ['TELESCOPE', 'DUT']]
        return [Plane(i, self.Config('TELESCOPE' if i < n_tel else 'DUT')) for i in range(n_tel + n_dut)]

    def init_residuals(self):
        from mod.residuals import ResidualAnalysis
        return ResidualAnalysis(self)

    def init_ref(self):
        from mod.reference import RefAnalysis
        return RefAnalysis(self)

    def init_tracks(self):
        from mod.tracks import TrackAnalysis
        return TrackAnalysis(self)

    def init_eff(self):
        from mod.efficiency import Efficiency
        return Efficiency(self)

    def init_tel(self):
        from mod.telescope import TelescopeAnalysis
        return TelescopeAnalysis(self)

    def load_file(self, test=False):
        if not test:
            self.Converter.run()
            try:
                f = h5py.File(self.Run.FileName, 'r')
                _ = f['Tracks']  # check if data is complete
                _ = f[str(self.Plane)]  # check if data is complete
                return f
            except (KeyError, OSError):
                warning(f'could not load data file {self.Run.FileName} -> start with dummy')
        return self.Dummy.load_file()

    def reload_data(self):
        self.F = self.load_file()

    def show_structure(self):
        print('.')
        for i0, (key, grp) in enumerate(self.F.items()):
            if 'Plane' not in key or key == 'Plane0' or key == f'Plane{self.Plane.Number}':
                print(f'├── {key}')
                for k, g in grp.items():
                    print(f'{"│" if i0 < len(self.F.keys()) - 1 else " "}   ├── {k}')
                    if hasattr(g, 'keys'):
                        for i in g.keys():
                            print(f'│   │   ├── {i}')
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA
    def get_time(self, cut=None):
        return self.get_data('Time', cut=cut, main_grp='Event').astype('f8')

    def get_data(self, grp, key=None, cut=None, pl=None, main_grp=None):
        data = self.F[choose(main_grp, str(self.Planes[choose(pl, self.Plane.Number)]))][grp]
        data = array(data) if key is None else array(data[key])
        return data if type(cut) is bool else self.Cut(cut, data, pl)

    def get_phs(self, e=False, cut=None):
        return self.get_data('Clusters', 'Charge', cut) * (self.DUT.VcalToEl if e else 1)

    def get_ph(self, cut=None):
        values = self.get_phs(cut=cut)
        m, s = mean_sigma(values)
        return ufloat(m, s / sqrt(values.size))

    def get_x(self, cut=None, pl=None):
        return self.get_data('Clusters', 'X', cut, pl)

    def get_y(self, cut=None, pl=None):
        return self.get_data('Clusters', 'Y', cut, pl)

    def get_tx(self, cut=None, pl=None):
        return self.get_data('Tracks', 'X', cut, pl)

    def get_ty(self, cut=None, pl=None):
        return self.get_data('Tracks', 'Y', cut, pl)

    def get_u(self, cut=None, pl=None, centre=False):
        return self.get_uv(cut, pl, centre)[0] if self.T else self.get_data('Clusters', 'U', cut, pl)

    def get_v(self, cut=None, pl=None, centre=False):
        return self.get_uv(cut, pl, centre)[1] if self.T else self.get_data('Clusters', 'V', cut, pl)

    def get_tu(self, cut=None, pl=None):
        return self.get_tuv(cut, pl)[0] if self.T else self.get_data('Tracks', 'U', cut, pl)

    def get_tv(self, cut=None, pl=None):
        return self.get_tuv(cut, pl)[1] if self.T else self.get_data('Tracks', 'V', cut, pl)

    def get_xy(self, local=True, cut=None, pl=None, centre=False):
        return array([self.get_x(cut, pl), self.get_y(cut, pl)]) if local else self.get_uv(cut, pl, centre)

    def get_uv(self, cut=None, pl=None, centre=False):
        return self.l2g(self.get_x(cut, pl), self.get_y(cut, pl), pl, centre) if self.T else array([self.get_u(cut, pl), self.get_v(cut, pl)])

    def get_txy(self, local=True, cut=None, pl=None, centre=False):
        return array([self.get_tx(cut, pl), self.get_ty(cut, pl)]) if local else self.get_tuv(cut, pl, centre)

    def get_tuv(self, cut=None, pl=None, centre=False):
        return self.l2g(self.get_tx(cut, pl), self.get_ty(cut, pl), pl, centre) if self.T else array([self.get_tu(cut, pl), self.get_tv(cut, pl)])

    def get_mask(self):
        return self.get_data('Mask', cut=False)

    def get_trigger_phase(self, cut=None, pl=None):
        return self.get_data('TriggerPhase', cut=cut, pl=pl)

    def get_cluster_size(self, cut=None, pl=None):
        return self.get_data('Clusters', 'Size', cut, pl)

    def get_chi2(self, cut=None):
        return self.get_data('Chi2', cut=cut, main_grp='Tracks') / self.get_data('Dof', cut=cut, main_grp='Tracks')
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region MISC
    def plane(self, pl):
        return self.Planes[choose(pl, self.Plane.Number)]

    @property
    def n(self):
        return self.get_x(cut=False).size

    def activate_surface(self, on=True):
        self.Cut.register('fid', self.Cut.make_fiducial(surface=on), 10, 'fid cut')
        self.Tracks.Cut.register('fid', self.Tracks.Cut.make_fiducial(surface=on), 30)
        self.Surface = on

    def deactivate_surface(self):
        self.activate_surface(False)

    def get_start_time(self):
        return datetime.fromtimestamp(self.Run.StartTime)

    def get_end_time(self):
        return datetime.fromtimestamp(self.Run.EndTime)

    def t_args(self, rel_t=False):
        return {'x_tit': 'Time [hh:mm]', 't_ax_off': self.Run.StartTime if rel_t else 0}

    @staticmethod
    def ax_tits(local=True):
        return {'x_tit': 'Column', 'y_tit': 'Row'} if local else {'x_tit': 'X [mm]', 'y_tit': 'Y [mm]'}

    @property
    def loc_bins(self):
        return bins.get_local(self.Plane)

    @property
    def ph_tit(self):
        return 'Pulse Height [vcal]'

    def get_ph_tit(self, e=False):
        return 'Charge [e]' if e else self.ph_tit

    def get_segments(self, nx, ny, width=False):
        x0, x1, y0, y1 = self.Cut.get_config('full size')
        if width:
            return arange(x0, x1 + (x1 - x0) // nx, nx, dtype='u2'), arange(y0, y1 + (y1 - y0) // ny, ny, dtype='u2')
        return linspace(x0, x1, nx + 1), linspace(y0, y1, ny + 1)

    def get_segment_centres(self, nx, ny, width=False):
        x, y = self.get_segments(nx, ny, width)
        return x[:-1] + diff(x) / 2, y[:-1] + diff(y) / 2

    def expand_inpixel(self, x, y, e=None, cell=False):
        cx, cy = [self.DUT.CellSize / 1000. / self.Plane.PX, self.DUT.CellSize / 1000. / self.Plane.PY]
        x, y = [(x % cx) / cx,  (y % cy) / cy] if cell else [x % 1, y % 1]
        # add edges of the neighbouring pixels
        xf = concatenate([x, x[x < .5] + 1, x[x >= .5] - 1, x[y < .5], x[y > .5]])
        yf = concatenate([y, y[x < .5], y[x >= .5], y[y < .5] + 1, y[y >= .5] - 1])
        ef = concatenate([e, e[x < .5], e[x >= .5], e[y < .5], e[y >= .5]]) if e is not None else None
        # add corners
        xf = concatenate([xf, x[(x < .5) & (y < .5)] + 1, x[(x >= .5) & (y >= .5)] - 1, x[(x < .5) & (y >= .5)] + 1, x[(x >= .5) & (y < .5)] - 1])
        yf = concatenate([yf, y[(x < .5) & (y < .5)] + 1, y[(x >= .5) & (y >= .5)] - 1, y[(x < .5) & (y >= .5)] - 1, y[(x >= .5) & (y < .5)] + 1])
        ef = concatenate([ef, e[(x < .5) & (y < .5)], e[(x >= .5) & (y >= .5)], e[(x < .5) & (y >= .5)], e[(x >= .5) & (y < .5)]]) if e is not None else None
        return [xf, yf] if e is None else [xf, yf, ef]
    # endregion MISC
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_mask(self, **dkw):
        masks = self.Cut.get_config('mask', default=zeros((2, 0))), self.get_mask(), self.Cut.get_thresh_mask(), self.Cut.get_cal_chi2_mask()
        d = concatenate([tile(m, i) for i, m in enumerate(masks, 1)], axis=1)
        self.Draw.histo_2d(*d, self.loc_bins, 'PixMask', **prep_kw(dkw, **self.ax_tits(), leg=self.Cut.get_fid(), pal=set_n_palette(len(masks)), rm=.03, stats=False, draw_opt='col', z_range=[0, 5]))
        Draw.legend([Draw.box(*[-1] * 4, line_color=0, fillcolor=i, show=False) for i in Draw.get_colors(len(masks))], ['Config', 'Noise Scan', 'High Thresh', 'Bad Calibration'], 'f', bottom=True)

    def draw_occupancy(self, local=True, bw=1, cut=None, fid=False, pl=None, **dkw):
        x, y = self.get_xy(local, self.Cut.get_nofid(cut, fid), pl)
        pl = self.Plane if local else self.Planes[0]
        return self.Draw.histo_2d(x, y, bins.get_xy(local, pl, bw, aspect_ratio=True), 'ClusterOcc', **prep_kw(dkw, qz=.99, z0=0, **self.ax_tits(local)))

    def draw_hit_map(self, bw=.3, local=True, cut=False, fid=False, **dkw):
        self.Tracks.draw_map(bw, local, self.Cut.get_nofid(cut, fid), local, **prep_kw(dkw, leg=self.Cut.get_fid() if local else None, title='HitMap'))

    def draw_cluster_size(self, cut=None, pl=None, **dkw):
        v = self.get_cluster_size(cut, pl)
        self.Draw.distribution(v, **prep_kw(dkw, title='Cluster Size', w=1, x0=-.5, q=1e-3, x_tit='Cluster Size'))

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        cut = self.Cut.get_nofid(cut, fid)
        (x, y), cs = self.Tracks.get_xy(local, cut), self.get_cluster_size(cut)
        self.Draw.prof2d(x, y, cs, bins.get_xy(local, self.Plane, res), 'Cluster Size', **prep_kw(dkw, z_tit='Cluster Size', **self.ax_tits(local)))

    def draw_trigger_phase(self, cut=None, **dkw):
        cut = -1 if type(cut) is bool else self.Cut.exclude('tp', cut)
        h = self.Draw.distribution(self.get_trigger_phase(cut), bins.TP, **prep_kw(dkw, title='Trigger Phase', x_tit='Trigger Phase'))
        format_histo(h, y_range=[0, 1.1 * h.GetMaximum()])

    def draw_time(self, **dkw):
        t = self.get_time()
        g = self.Draw.profile(arange(t.size), t, bins.make(0, t.size, sqrt(t.size)), 'Time', **prep_kw(dkw, markersize=.6, x_tit='Event Number', y_tit='Time [hh:mm]', draw_opt='aplx', graph=True))
        set_time_axis(g, axis='Y')

    def draw_grid(self, nx=2, ny=3, w=1, width=False):
        self.Draw.grid(*self.get_segments(nx, ny, width), w)

    def draw_inpixel_map(self, res=.1, cut=None, cell=False, show=True):
        x, y = self.expand_inpixel(cell=cell, *self.Tracks.get_xy(cut=cut))
        self.Draw.histo_2d(x, y, bins.get_pixel(self.Plane, res, cell=cell), 'Hit Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region CORRELATION
    def draw_x_correlation(self, pl=2, pl1=None, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        return self.Draw.histo_2d(self.get_x(c, pl), self.get_x(c, pl1), **prep_kw(dkw, title='XCorr', x_tit=f'Column Plane {pl}', y_tit=f'Column Plane {choose(pl1, self.Plane.Number)}'))

    def draw_y_correlation(self, pl=2, pl1=None, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        return self.Draw.histo_2d(self.get_y(c, pl), self.get_y(c, pl1), **prep_kw(dkw, title='YCorr', x_tit=f'Row Plane {pl}', y_tit=f'Row Plane {choose(pl1, self.Plane.Number)}'))

    def draw_correlation_trend(self, pl=0, pl1=None, thresh=.2, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        d0, d1, t = self.get_xy(pl=pl, cut=c), self.get_xy(pl=pl1, cut=c), self.get_time(c)
        g = [self.Draw.graph(*get_3d_correlations(self.Draw.histo_3d(t, d0[i], d1[i]), thresh=thresh), y_tit='Correlation Factor', show=False) for i in range(2)]
        return self.Draw.multigraph(g, 'CorrFac', ['x', 'y'], draw_opt='pl', **prep_kw(dkw, **self.t_args(), y_range=[-1.05, 1.05]))

    def draw_alignment(self, pl=2, thresh=.3, **dkw):
        gx, gy = self.draw_correlation_trend(pl, show=False).GetListOfGraphs()
        (t, x), y = get_graph_vecs(gx, err=False), get_graph_y(gy, err=False)
        r = [1 if abs(ix) > thresh and abs(iy) > thresh else 2 for ix, iy in zip(x, y)]
        x, y = t.repeat(r), ones(sum(r))
        binning = bins.from_vec(get_graph_x(gx)) + [3, 0, 3]
        gStyle.SetPalette(3, array([1, 633, 418], 'i'))
        self.Draw.histo_2d(x, y, binning, 'Event Alignment', **prep_kw(dkw, **self.t_args(), y_tit='Alignment', stats=False, l_off_y=99, center_y=True, draw_opt='col', z_range=[0, 2]))
        Draw.legend([Draw.box(0, 0, 0, 0, line_color=c, fillcolor=c) for c in [418, 633]], ['aligned', 'misaligned'], 'f')
    # endregion CORRELATION
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL
    def draw_signal_distribution(self, cut=None, draw_thresh=False, e=False, **dkw):
        trim = self.Calibration.Trim * (self.DUT.VcalToEl if e else 1)
        legs = [self.Draw.vertical_line(trim, 0, 1e5, w=2, show=False), self.Draw.tlatex(.2, .5, f'Threshold = {trim:.0f} vcal', ndc=True, angle=90, size=.04, show=False)] if draw_thresh else None
        return self.Draw.distribution(self.get_phs(e, cut), **prep_kw(dkw, title='PH', x_tit=self.ph_tit, leg=legs))

    def draw_charge_distribution(self, cut=None, draw_thresh=False, **dkw):
        return self.draw_signal_distribution(cut, draw_thresh, e=True, **dkw)

    def draw_low_ph_map(self, cmax, cmin=None, res=.5, **dkw):
        self.draw_hit_map(res, cut=self.Cut.get_nofid() + self.Cut.make_ph(cmax, cmin), **dkw)

    def draw_signal_map(self, res=.3, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [self.Tracks.get_xy, self.get_phs]]
        self.Draw.prof2d(x, y, z_, bins.get_local(self.Plane, res), 'Charge Map', **prep_kw(dkw, leg=self.Cut.get_fid(), z_tit=self.ph_tit, **self.ax_tits()))

    def draw_signal_occupancy(self, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [self.get_xy, self.get_phs]]
        self.Draw.prof2d(x, y, z_, self.loc_bins, 'Charge Occupancy', **prep_kw(dkw, leg=self.Cut.get_fid(), z_tit=self.ph_tit, **self.ax_tits()))

    def draw_signal_vs_tp(self, cut=None, **dkw):
        x, y = [f(cut=self.Cut.exclude('tp', cut)) for f in [self.get_trigger_phase, self.get_phs]]
        self.Draw.profile(x, y, bins.TP, 'Charge vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', y_tit=self.ph_tit, graph=True))

    def draw_signal_trend(self, bw=None, e=False, cut=None, **dkw):
        x, y = self.get_time(cut), self.get_phs(e=e, cut=cut)
        g = self.Draw.profile(x, y, find_bins(x, w=bw), **self.t_args(), graph=True, y_tit=self.get_ph_tit(e), stats=True, show=False)
        return self.Draw(g, **prep_kw(dkw, y_range=ax_range(get_graph_y(g, err=False), fl=1, fh=2)))

    def fit_signal(self, bw=None, e=False, **dkw):
        g = self.draw_signal_trend(bw, e, show=False)
        fit = FitRes(g.Fit('pol0', 'sq'))
        self.Draw(g, **prep_kw(dkw, stats=set_statbox(fit=True), show=False))
        return fit

    def draw_inpixel_charge(self, res=.1, cut=None, show=True, cell=False):
        (x, y), c = self.Tracks.get_xy(cut=cut), self.get_phs(cut=cut)
        x, y, c = self.expand_inpixel(x, y, c, cell)
        self.Draw.prof2d(x, y, c, bins.get_pixel(self.Plane, res, cell=cell), 'Charge Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
    # endregion SIGNAL
    # ----------------------------------------

    # ----------------------------------------
    # region COORDINATE TRANSFORM
    @property
    def inv_x(self):
        return self.Converter.get_alignment()['sensors'][self.Plane.Number]['unit_u'][0] < 0

    def l2g(self, x, y, pl=None, centre=False, inv_x=None, invert=False):
        pl = self.plane(pl)
        a = self.Converter.get_alignment()['sensors'][pl.Number]
        ox, oy = array(a['offset'][:2]) - (array([pl.W, pl.H]) / 2 if centre else 0)
        rx, ry = array(a['unit_u']) * (-1 if choose(inv_x, self.inv_x) else 1), a['unit_v']
        return transform(x, y, sx=pl.PX, sy=pl.PY, ox=ox, oy=oy, rx=rx, ry=ry, order='trs', invert=invert)

    def g2l(self, x, y, pl=None, centre=False, inv_x=None):
        return self.l2g(x, y, pl, centre, inv_x, invert=True)
    # endregion COORDINATE TRANSFORM
    # ----------------------------------------

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(show=show)
        get_last_canvas().Modified()
        get_last_canvas().Update()
        if fit.get_chi2() > chi_thresh and nconv < 80:
            Draw.Count += 5
            self.info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.get_chi2()))
            fit = self.fit_langau(h, nconv + Draw.Count, chi_thresh=chi_thresh, show=show)
        print('MPV: {:1.1f}'.format(fit.get_mpv()))
        Draw.Count = 0
        self.Draw.add(fit)
        return fit
