#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a single DUT
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
import cern.converter
import src.bins as bins
import src.converter
from mod.dut_cuts import DUTCut
from plotting.fit import *
from src.analysis import *
from src.currents import Currents
from src.run import Run
from src.dut import Plane
from utility.utils import *
from utility.affine_transformations import transform, m_transform


def no_trans(f):
    def inner(*args, **kwargs):
        old = DUTAnalysis.Trans
        DUTAnalysis.Trans = False
        v = f(*args, **kwargs)
        DUTAnalysis.Trans = old
        return v
    return inner


class DUTAnalysis(Analysis):

    Trans = True  # use internal algorithm to improve alignment of the local track coordinates
    L2G = False   # transform local to global coordinates instead of using global directly

    def __init__(self, run_number, dut_number, test_campaign, verbose=True, test=False):

        Analysis.__init__(self, test_campaign, meta_sub_dir='DUT', verbose=verbose)
        self.Run = Run.from_ana(run_number, dut_number, self)
        self.print_start()

        # data
        self.DUT = self.Run.DUT
        self.Converter = self.init_converter()
        self.Proteus = self.Converter.Proteus
        self.Planes = self.init_planes()
        self.Plane = self.Planes[self.DUT.Plane.Number]  # update rotated

        if test:
            return

        self.F = self.load_file()

        # INFO
        self.N = self.n
        self.NEvents = self.F['Event']['Time'].size
        self.NTracks = self.F['Tracks']['Size'].size
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()
        self.Duration = (self.EndTime - self.StartTime).seconds
        self.Surface = False

        # SUBCLASSES
        self.Calibration = self.Converter.load_calibration(self.DUT.Number)
        self.Cut = DUTCut(self)
        self.Residuals = self.init_residuals()

        self.Tel = self.init_tel()
        self.REF = self.init_ref()
        self.Tracks = self.init_tracks()

        self.Cut.make_additional()
        self.Efficiency = self.init_eff()
        self.Resolution = self.init_resolution()
        self.Currents = Currents(self)

        if self.BeamTest.Location == 'CERN':
            self.verify_alignment()

    def __repr__(self):
        return f'{self} of {self.Run.DUT}, run {self.Run} ({self.BeamTest}), {self.ev_str}'

    @property
    def ev_str(self):
        return f'{ev2str(self.NEvents if hasattr(self, "NEvents") else self.Run.n_ev)} ev'

    @classmethod
    def from_run(cls, run: Run, verbose=True, test=False):
        return cls(run.Number, run.DUT.Number, run.TCDir.stem, verbose, test)

    @property
    def suffix(self):
        return f'{self.DUT}-{self.Run}-{self.BeamTest.Location}'.lower().replace('ii6-', '')

    def make_plots(self, res=.2, n=50):
        old_dir = self.Draw.ResultsDir
        SaveDraw.SaveOnServer = False
        self.Draw.ResultsDir = Dir.joinpath('tmp')
        rx, ry = self.Cut.get_config('dia size').reshape(2, -1) if 'dia size' in self.Cut.Config.options() else (None, None)
        self.draw_signal_distribution(x0=-50, x1=1040, fn=f'sd-{self.suffix}')
        for pal in [53, 55]:
            self.draw_occupancy(cut=0, x_range=rx, y_range=ry, pal=pal, fn=f'occ-{self.suffix}-p{pal}')
            self.Efficiency.draw_map(res=res, x_range=rx, y_range=ry, pal=pal, fn=f'em-{self.suffix}-p{pal}', leg=self.Cut.get_fid())
            self.draw_signal_map(res=res, x_range=rx, y_range=ry, pal=pal, qz=.999, fn=f'sm-{self.suffix}-p{pal}')
            self.draw_cluster_size_map(res=.1, x_range=rx, y_range=ry, pal=pal, qz=.995, fn=f'csm-{self.suffix}-p{pal}')
            self.draw_ph_in_pixel(n, 75, 50, pal=pal, fn=f'ph-pix-{self.suffix}-p{pal}')
            self.draw_cs_in_pixel(n, 75, 50, pal=pal, qz=.999, fn=f'cs-pix-{self.suffix}-p{pal}')
            self.Efficiency.draw_in_pixel(n, 75, 50, pal=pal, z_range=[70, 100], fn=f'e-pix-{self.suffix}-p{pal}')
            self.draw_hitmap_in_pixel(n, 75, 50, pal=pal, fn=f'hm-pix-{self.suffix}-p{pal}')
        self.Draw.ResultsDir = old_dir
        SaveDraw.SaveOnServer = True

    # ----------------------------------------
    # region INIT
    @property
    def server_save_dir(self):
        return Path('duts', str(self.DUT), self.BeamTest.Tag, str(self.Run))

    def init_converter(self):
        return self.converter.from_run(self.Run)

    @property
    def converter(self):
        return cern.converter.CERNConverter if self.BeamTest.Location == 'CERN' else src.converter.Converter

    def init_planes(self):
        n_tel, n_dut = self.Converter.NTelPlanes, self.Converter.NDUTPlanes
        rot = [abs(s['unit_u'][1]) > .5 for s in self.Proteus.alignment()['sensors']] if self.Proteus.has_alignment else [False] * (n_tel + n_dut + self.Proteus.NRefPlanes)
        if len(rot) != n_tel + n_dut + self.Proteus.NRefPlanes:
            warning(f'Number of sensors in alignment ({len(rot)}) does not agree with found sensors ({n_tel + n_dut + ("REF" in Analysis.Config)})')
            return [Plane(i) for i in range(10)]
        pl = [Plane(i, typ='TELESCOPE', rotated=rot[i]) for i in range(n_tel)]
        pl += [Plane(n_tel, typ='REF', rotated=rot[n_tel])] if self.Proteus.NRefPlanes > 0 else []
        return pl + [Plane(len(pl) + i, typ='DUT', rotated=rot[len(pl) + i]) for i in range(n_dut)]

    def init_residuals(self):
        from mod.residuals import res_analysis
        return res_analysis(self.__class__)(self)

    def init_ref(self):
        if self.Run.NDUTs > 1 or self.Proteus.NRefPlanes:
            from mod.reference import ref_analysis
            return ref_analysis(self.__class__)(self)

    def init_tracks(self):
        from mod.tracks import track_analysis
        return track_analysis(self.__class__)(self)

    def init_eff(self):
        from mod.efficiency import eff_analysis
        return eff_analysis(self.__class__)(self)

    def init_tel(self):
        from mod.telescope import tel_analysis
        return tel_analysis(self.__class__)(self)

    def init_resolution(self):
        if self.REF is not None:
            from mod.resolution import reso_analysis
            return reso_analysis(self.__class__)(self.REF)

    @property
    def file_name(self):
        return self.Run.FileName

    def remove_file(self):
        remove_file(self.file_name)

    def load_file(self):
        if self.Converter.run():
            self.remove_metadata()
        try:
            f = h5py.File(self.file_name, 'r')
            _ = f['Tracks']         # check if data is complete
            _ = f[str(self.Plane)]  # check if data is complete
            return f
        except (KeyError, OSError) as err:
            self.remove_file()
            critical(f'error loading data file, deleting {self.file_name}\n{err}')

    def reload_data(self):
        self.F = self.load_file()

    def show_structure(self):
        show_hdf5(self.F, str(self.Plane), 'Plane0', ex_str='Plane')

    def has_alignment(self, imax=20):
        return all([file_exists(self.make_pickle_path('AM', imax, 'alignment', run='', dut=i)) for i in [self.DUT.Number, self.REF.DUT.Number]])

    def verify_alignment(self, cut=100):
        """checks if the residuals are centred around 0 for a given cut [um]"""
        mx, my = self.Residuals.means() * 1e3  # mm -> um
        if not sqrt(mx.n ** 2 + my.n ** 2) < cut:
            self.Converter.realign()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA
    def time(self, cut=None):
        return self.get_data('Time', cut=cut, main_grp='Event').astype('f8')

    def mean_time(self):
        t0, t1 = self.F['Event']['Time'][[0, -1]]
        return ufloat(mean([t0, t1]), (t1 - t0) / 2)

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

    def get_x(self, cut=None, pl=None, rot=False):
        return self.get_y(cut, pl) if self.plane(pl).Rotated and rot else self.get_data('Clusters', 'X', cut, pl)

    def get_y(self, cut=None, pl=None, rot=False):
        return self.get_x(cut, pl) if self.plane(pl).Rotated and rot else self.get_data('Clusters', 'Y', cut, pl)

    def get_tx(self, cut=None, pl=None):
        return self.get_data('Tracks', 'X', cut, pl)

    def get_ty(self, cut=None, pl=None):
        return self.get_data('Tracks', 'Y', cut, pl)

    def get_u(self, cut=None, pl=None, centre=False):
        return self.get_uv(cut, pl, centre)[0] if DUTAnalysis.L2G else self.get_data('Clusters', 'U', cut, pl)

    def get_v(self, cut=None, pl=None, centre=False):
        return self.get_uv(cut, pl, centre)[1] if DUTAnalysis.L2G else self.get_data('Clusters', 'V', cut, pl)

    def get_tu(self, cut=None, pl=None):
        return self.get_tuv(cut, pl)[0] if DUTAnalysis.L2G else self.get_data('Tracks', 'U', cut, pl)

    def get_tv(self, cut=None, pl=None):
        return self.get_tuv(cut, pl)[1] if DUTAnalysis.L2G else self.get_data('Tracks', 'V', cut, pl)

    def get_xy(self, local=True, cut=None, pl=None, centre=False, rot=False):
        return array([self.get_x(cut, pl, rot), self.get_y(cut, pl, rot)]) if local else self.get_uv(cut, pl, centre)

    def get_uv(self, cut=None, pl=None, centre=False):
        return self.l2g(self.get_x(cut, pl), self.get_y(cut, pl), pl, centre) if DUTAnalysis.L2G else array([self.get_u(cut, pl), self.get_v(cut, pl)])

    def get_txy(self, local=True, cut=None, pl=None, centre=False):
        d = array([self.get_tx(cut, pl), self.get_ty(cut, pl)]) if local else self.get_tuv(cut, pl, centre)
        return m_transform(self.Residuals.m, *d) if DUTAnalysis.Trans and local else d

    def get_tuv(self, cut=None, pl=None, centre=False):
        return self.l2g(self.get_tx(cut, pl), self.get_ty(cut, pl), pl, centre) if DUTAnalysis.L2G else array([self.get_tu(cut, pl), self.get_tv(cut, pl)])

    def get_mask(self):
        return self.get_data('Mask', cut=False)

    def get_trigger_phase(self, cut=None, pl=None):
        return self.get_data('Trigger', 'Phase', cut=cut, pl=pl)

    def get_cluster_size(self, cut=None, pl=None):
        return self.get_data('Clusters', 'Size', cut, pl)

    def get_chi2(self, cut=None):
        return self.get_data('Chi2', cut=cut, main_grp='Tracks') / self.get_data('Dof', cut=cut, main_grp='Tracks')

    def get_slope_x(self, cut=None):
        return self.get_data('SlopeX', cut=cut, main_grp='Tracks')

    def get_slope_y(self, cut=None):
        return self.get_data('SlopeY', cut=cut, main_grp='Tracks')

    def alignment(self, pl):
        return self.Proteus.alignment()['sensors'][pl]

    def eff(self):
        return self.Efficiency.value()

    def current(self):
        return self.Currents()
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

    def segments(self, nx, ny, width=False):
        x0, x1, y0, y1 = self.Cut.get_config('full size')
        if width:
            return arange(x0, x1 + (x1 - x0) // nx, nx, dtype='u2'), arange(y0, y1 + (y1 - y0) // ny, ny, dtype='u2')
        return linspace(x0, x1, nx + 1), linspace(y0, y1, ny + 1)

    def segment_centres(self, nx, ny, width=False):
        x, y = self.segments(nx, ny, width)
        return x[:-1] + diff(x) / 2, y[:-1] + diff(y) / 2
    # endregion MISC
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_mask(self, **dkw):
        masks = self.Cut.get_config('mask', default=zeros((2, 0))), self.get_mask(), self.Cut.get_thresh_mask(), self.Cut.get_cal_chi2_mask()
        d = concatenate([tile(m, i) for i, m in enumerate(masks, 1)], axis=1)
        self.Draw.histo_2d(*d, self.loc_bins, 'PixMask', **prep_kw(dkw, **self.ax_tits(), leg=self.Cut.get_fid(), pal=n_pal(len(masks)), rm=.03, stats=False, draw_opt='col', z_range=[0, 5]))
        Draw.legend([Draw.box(*[-1] * 4, line_color=0, fillcolor=i, show=False) for i in Draw.get_colors(len(masks))], ['Config', 'Noise Scan', 'High Thresh', 'Bad Calibration'], 'f', bottom=True)
        self.Draw.save_plots('Mask')

    def draw_occupancy(self, local=True, bw=1, cut=None, fid=False, pl=None, **dkw):
        x, y = self.get_xy(local, self.Cut.get_nofid(cut, fid), pl)
        pl = self.Plane if local else self.Planes[0]
        return self.Draw.histo_2d(x, y, bins.get_xy(local, pl, bw, aspect_ratio=True), 'ClusterOcc', **prep_kw(dkw, qz=.99, z0=0, **self.ax_tits(local), file_name='Occupancy'))

    def draw_hit_map(self, bw=.3, local=True, cut=False, fid=False, **dkw):
        return self.Tracks.draw_map(bw, local, self.Cut.get_nofid(cut, fid), local, **prep_kw(dkw, leg=self.Cut.get_fid() if local else None, title='HitMap', file_name='HitMap'))

    def draw_cluster_size(self, cut=None, pl=None, **dkw):
        v = self.get_cluster_size(self.Cut.exclude('cs', cut), pl)
        self.Draw.distribution(v, **prep_kw(dkw, title='Cluster Size', w=1, x0=-.5, q=1e-3, x_tit='Cluster Size', file_name='ClusterSize'))

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        cut = self.Cut.get_nofid(self.Cut.exclude('cs', cut), fid)
        (x, y), cs = self.get_txy(local, cut), self.get_cluster_size(cut)
        self.Draw.prof2d(x, y, cs, bins.get_xy(local, self.Plane, res), 'Cluster Size', **prep_kw(dkw, qz=.98, z0=1, z_tit='Cluster Size', **self.ax_tits(local), file_name='CSMap'))

    def draw_trigger_phase(self, cut=None, **dkw):
        h = self.Draw.distribution(self.get_trigger_phase(self.Cut.exclude('tp', cut)), bins.TP, **prep_kw(dkw, title='Trigger Phase', x_tit='Trigger Phase'))
        self.Draw(h, **prep_kw(dkw, y_range=[0, 1.1 * h.GetMaximum()], file_name='TriggerPhase'))

    def draw_time(self, cut=None, **dkw):
        t = self.time(cut)
        g = self.Draw.profile(arange(t.size), t, title='Time', **prep_kw(dkw, markersize=.6, x_tit='Event Number', y_tit='Time [hh:mm]', draw_opt='aplx', graph=True))
        set_time_axis(g, axis='Y')
        self.Draw.save_plots('Time')

    def draw_time_dist(self, cut=None, **dkw):
        self.Draw.distribution(self.time(cut), **prep_kw(dkw, title='TimeDist', **self.t_args(), stats=set_statbox(entries=True)))

    def draw_grid(self, nx=2, ny=3, w=1, width=False):
        return self.Draw.grid(*self.segments(nx, ny, width), width=w)

    def draw_geo(self, top=True):
        x, y = self.Proteus.z_positions(raw=True), [pl.W if top else pl.H for pl in self.Planes]
        planes = [self.Draw.vertical_line(x, *array([-.5, .5]) * y, show=False) for x, y in zip(x, y)]
        tits = {'y_tit': 'y [mm]' if top else 'x [mm]', 'file_name': f'Geo{"Top" if top else "Side"}'}
        return self.Draw.graph([-.05 * max(x), max(x)], array([-.6, .6]) * max(y), x_tit='z [mm]', **tits, markersize=0, leg=planes)

    def draw_geo_side(self):
        return self.draw_geo(top=False)
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region CORRELATION
    def draw_u_corr(self, pl=2, pl1=None, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        x = [self.get_u(c, p) for p in [pl, pl1]]
        return self.Draw.histo_2d(*x, **prep_kw(dkw, title='UCorr', x_tit=f'X P{pl} [#mum]', y_tit=f'Y P{self.plane(pl1).Number} [#mum]', file_name='UCorr'))

    def draw_x_correlation(self, pl=2, pl1=None, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        x = [self.get_x(c, p, rot=True) for p in [pl, pl1]]
        return self.Draw.histo_2d(*x, **prep_kw(dkw, title='XCorr', x_tit=f'Column Plane {pl}', y_tit=f'Column Plane {choose(pl1, self.Plane.Number)}', file_name='XCorr'))

    def draw_y_correlation(self, pl=2, pl1=None, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        y = [self.get_y(c, p, rot=True) for p in [pl, pl1]]
        return self.Draw.histo_2d(*y, **prep_kw(dkw, title='YCorr', x_tit=f'Row Plane {pl}', y_tit=f'Row Plane {choose(pl1, self.Plane.Number)}', file_name='YCorr'))

    def draw_correlation_trend(self, pl=0, pl1=None, thresh=.2, **dkw):
        c = self.Cut.make_correlation(pl, pl1)
        d0, d1, t = self.get_xy(pl=pl, cut=c, rot=True), self.get_xy(pl=pl1, cut=c, rot=True), self.time(c)
        g = [self.Draw.graph(*get_3d_correlations(self.Draw.histo_3d(t, d0[i], d1[i]), thresh=thresh), y_tit='Correlation Factor', show=False) for i in range(2)]
        return self.Draw.multigraph(g, 'CorrFac', ['x', 'y'], draw_opt='pl', **prep_kw(dkw, **self.t_args(), y_range=[-1.05, 1.05], file_name='CorrTrend'))

    def draw_alignment(self, pl=2, thresh=.3, **dkw):
        gx, gy = self.draw_correlation_trend(pl, show=False).GetListOfGraphs()
        (t, x), y = get_graph_vecs(gx, err=False), get_graph_y(gy, err=False)
        r = [1 if abs(ix) > thresh and abs(iy) > thresh else 2 for ix, iy in zip(x, y)]
        x, y = t.repeat(r), ones(sum(r))
        binning = bins.from_vec(get_graph_x(gx)) + [3, 0, 3]
        gStyle.SetPalette(3, array([1, 633, 418], 'i'))
        self.Draw.histo_2d(x, y, binning, 'Event Alignment', **prep_kw(dkw, **self.t_args(), y_tit='Alignment', stats=False, l_off_y=99, center_y=True, draw_opt='col', z_range=[0, 2]))
        Draw.legend([Draw.box(0, 0, 0, 0, line_color=c, fillcolor=c) for c in [418, 633]], ['aligned', 'misaligned'], 'f')
        self.Draw.save_plots('EventAlignment')
    # endregion CORRELATION
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL
    def draw_signal_distribution(self, cut=None, draw_thresh=False, e=False, **dkw):
        return self.Draw.distribution(self.get_phs(e, cut), **prep_kw(dkw, title='PH', x_tit=self.ph_tit, leg=self.draw_trim(e, draw_thresh), file_name='SignalDist'))

    def draw_trim(self, e, thresh=False):
        if self.Calibration.Trim is None:
            return
        trim = self.Calibration.Trim * (self.DUT.VcalToEl if e else 1)
        return [self.Draw.vertical_line(trim, 0, 1e5, w=2, show=False), self.Draw.tlatex(.2, .5, f'Threshold = {trim:.0f} vcal', ndc=True, angle=90, size=.04, show=False)] if thresh else None

    def draw_charge_distribution(self, cut=None, draw_thresh=False, **dkw):
        return self.draw_signal_distribution(cut, draw_thresh, e=True, **dkw)

    def draw_low_ph_map(self, cmax, cmin=None, res=.5, **dkw):
        self.draw_hit_map(res, cut=self.Cut.get_nofid() & self.Cut.make_ph(cmax, cmin), **prep_kw(dkw, file_name=f'LowPH{cmax}'))

    def draw_signal_map(self, res=.3, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [self.get_txy, self.get_phs]]
        return self.Draw.prof2d(x, y, z_, bins.get_local(self.Plane, res), 'Charge Map', **prep_kw(dkw, qz=.95, leg=self.Cut.get_fid(), z_tit=self.ph_tit, **self.ax_tits(), file_name='SignalMap'))

    def draw_signal_occupancy(self, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [self.get_xy, self.get_phs]]
        self.Draw.prof2d(x, y, z_, self.loc_bins, 'Charge Occupancy', **prep_kw(dkw, leg=self.Cut.get_fid(), z_tit=self.ph_tit, **self.ax_tits()))

    def draw_signal_vs_tp(self, cut=None, **dkw):
        x, y = [f(cut=self.Cut.exclude('tp', cut)) for f in [self.get_trigger_phase, self.get_phs]]
        self.Draw.profile(x, y, bins.TP, 'Charge vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', y_tit=self.ph_tit, graph=True, file_name='PH-TP'))

    def draw_signal_trend(self, bw=None, e=False, cut=None, **dkw):
        x, y = self.time(cut), self.get_phs(e=e, cut=cut)
        g = self.Draw.profile(x, y, find_bins(x, w=bw), **self.t_args(), graph=True, y_tit=self.get_ph_tit(e), stats=True, show=False)
        return self.Draw(g, **prep_kw(dkw, y_range=ax_range(get_graph_y(g, err=False), fl=1, fh=2)))

    def fit_signal(self, bw=None, e=False, **dkw):
        g = self.draw_signal_trend(bw, e, show=False)
        fit = FitRes(g.Fit('pol0', 'sq'))
        self.Draw(g, **prep_kw(dkw, stats=set_statbox(fit=True), show=False))
        return fit

    def draw_eff_vs_ph(self, res=.5):
        ph, e = self.draw_signal_map(res, fid=True), self.Efficiency.draw_map(res, fid=True)
        return self.Draw.maps_profile(ph, e, file_name='EffVsPh')
    # endregion SIGNAL
    # ----------------------------------------

    # ----------------------------------------
    # region COORDINATE TRANSFORM
    @property
    def inv_x(self):
        return self.Converter.Proteus.alignment()['sensors'][self.Plane.Number]['unit_u'][0] < 0

    def l2g(self, x, y, pl=None, centre=False, inv_x=None, invert=False):
        pl = self.plane(pl)
        a = self.Converter.Proteus.alignment()['sensors'][pl.Number]
        ox, oy = array(a['offset'][:2]) - (array([pl.W, pl.H]) / 2 if centre else 0)
        rx, ry = array(a['unit_u']) * (-1 if choose(inv_x, self.inv_x) else 1), a['unit_v']
        return transform(x, y, sx=pl.PX, sy=pl.PY, ox=ox, oy=oy, rx=rx, ry=ry, order='trs', invert=invert)

    def g2l(self, x, y, pl=None, centre=False, inv_x=None):
        return self.l2g(x, y, pl, centre, inv_x, invert=True)
    # endregion COORDINATE TRANSFORM
    # ----------------------------------------

    # ----------------------------------------
    # region IN PIXEL
    def contracted_vars(self, mx=1, my=1, ox=0, oy=0, fz=None, cut=None, contract=True):
        x, y = self.get_txy(cut=cut)
        z_ = self.get_phs(cut=cut) if fz is None else zeros(x.size, dtype='?') if not fz else fz(cut=cut)
        x, y = (x + ox / self.Plane.PX / 1e3) % mx, (y + oy / self.Plane.PY / 1e3) % my
        return array(self.expand_vars(x, y, z_, mx, my) if contract else (x, y, z_)) * [[self.Plane.PX * 1e3], [self.Plane.PY * 1e3], [1]]  # convert from pixel to um

    @staticmethod
    def expand_vars(x, y, z_, mx, my):
        """copy the vars x, y, z to half a cell of size [mx, my] in each direction"""
        d = array([x, y]).T
        (x, y), e = concatenate([d + [i, j] for i in [-mx, 0, mx] for j in [-my, 0, my]]).T, tile(z_, 9)  # copy arrays in each direction
        cut = (x >= -mx / 2) & (x <= mx * 3 / 2) & (y >= -my / 2) & (y <= my * 3 / 2)  # select only half of the copied cells
        return x[cut], y[cut], e[cut]

    def draw_in(self, mx, my, ox=0, oy=0, n=None, cut=None, fz=None, **dkw):
        x, y, z_ = self.contracted_vars(mx / self.Plane.PX * 1e-3, my / self.Plane.PY * 1e-3, ox, oy, fz, cut)
        n = choose(n, freedman_diaconis, x=x) // 2 * 2  # should be symmetric...
        d = lambda w: round((n + .5) * (max(mx, my) / n - w) / w) * w  # extra spacing to account for different mx and my
        binning = sum([make_bins(-(i + w) / 2 - d(w), (3 * i + w) / 2 + d(w), w, last=True) for i, w in [(mx, mx / n), (my, my / n)]], start=[])
        cell = self.Draw.box(0, 0, mx, my, width=2, show=False, fillstyle=1)
        fh = self.Draw.prof2d if any(z_) else self.Draw.histo_2d
        h = fh(x, y, zz=z_, binning=binning, save=False, show=False, **prep_kw(rm_key(dkw, 'show'), title='Signal In Cell', x_tit='X [#mum]', y_tit='Y [#mum]', z_tit='Pulse Height [vcal]'))
        return self.Draw(h, **prep_kw(dkw, leg=self.draw_columns(show=get_kw('show', dkw, default=True)) + [cell]))

    def draw_in_cell(self, ox=0, oy=0, n=None, cut=None, fz=None, tit='PH', **dkw):
        return self.draw_in(*self.DUT.PXY, ox, oy, n, cut, fz, **prep_kw(dkw, title=f'{tit} in Cell', file_name=f'{tit.title()}InCell'))

    def draw_in_pixel(self, ox=0, oy=0, n=None, cut=None, fz=None, tit='PH', **dkw):
        return self.draw_in(*self.Plane.PXY * 1e3, ox, oy, n, cut, fz, **prep_kw(dkw, title=f'{tit} in Pixel', file_name=f'{tit.title()}InPixel'))

    def draw_hitmap_in_pixel(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_pixel(ox, oy, n, cut, fz=False, tit='HitMap', **dkw)

    def draw_hitmap_in_cell(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_cell(ox, oy, n, cut, fz=False, tit='HitMap', **dkw)

    def draw_ph_in_cell(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_cell(ox, oy, n, cut, **prep_kw(dkw, pal=53))

    def draw_ph_in_pixel(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_pixel(ox, oy, n, cut, **prep_kw(dkw, pal=53))

    def draw_cs_in_cell(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_cell(ox, oy, n, self.Cut.exclude('cs', cut), fz=self.get_cluster_size, tit='CS', **prep_kw(dkw, qz=.98, z0=1, pal=53, z_tit='Cluster Size'))

    def draw_cs_in_pixel(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_pixel(ox, oy, n, cut, fz=self.get_cluster_size, tit='CS', **prep_kw(dkw, qz=.98, z0=1, pal=53, z_tit='Cluster Size'))

    def draw_tp_in_pixel(self, n=None, ox=0, oy=0, cut=None, **dkw):
        return self.draw_in_pixel(ox, oy, n, self.Cut.exclude('tp', cut), fz=self.get_trigger_phase, tit='TP', **prep_kw(dkw, pal=53, z_tit='Trigger Phase'))

    def draw_columns(self, show=True):
        if hasattr(self.DUT, 'ColumnDiameter'):
            wx, wy, c, d = self.DUT.PX, self.DUT.PY, get_last_canvas(), self.DUT.ColumnDiameter.n
            x0, x1, y0, y1 = c.GetUxmin(), c.GetUxmax(), c.GetUymin(), c.GetUymax()
            b = [Draw.circle(d / 2, x, y, fill_color=602, fill=True, show=show) for x in arange(-2 * wx, x1, wx) for y in arange(-2 * wy, y1, wy) if x > x0 and y > y0]      # bias
            r = [Draw.circle(d / 2, x, y, fill_color=799, fill=True, show=show) for x in arange(-2.5 * wx, x1, wx) for y in arange(-2.5 * wy, y1, wy) if x > x0 and y > y0]  # readout
            g = [Draw.make_tgraph([1e3], [1e3], color=i, show=False, markersize=2) for i in [602, 799]]  # dummy graphs for legend
            return [Draw.legend(g, ['bias', 'readout'], 'p', y2=.82, show=show)] + b + r
        return []
    # endregion IN PIXEL
    # ----------------------------------------

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(draw=show)
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
