#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a single DUT
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F
from numpy import invert, arctan

from src.analysis import *
from src.calibration import Calibration
from src.cern_run import CERNRun
from src.converter import Converter
from src.currents import Currents
from src.cut import Cut
from src.desy_converter import DESYConverter
from src.desy_run import DESYRun
from src.dummy import Dummy
from src.dut import Plane
from plotting.fit import *
from src.reference import RefAnalysis
from src.utils import *
from src.telescope import TelescopeAnalysis
from src.tracks import TrackAnalysis
import src.bins as bins


# TODO: add residuals to data!
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
        self.Data = self.load_data(test)
        self.init_cuts()

        # INFO
        self.NEvents = self.get_entries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()
        self.Duration = (self.EndTime - self.StartTime).seconds
        self.Surface = False

        # SUBCLASSES
        self.Calibration = Calibration(self.Run)
        self.Tel = TelescopeAnalysis(self)
        self.Tracks = TrackAnalysis(self)
        self.REF = RefAnalysis(self)
        self.Currents = Currents(self)
        self.add_cuts()
        self.add_track_cuts()

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

    def init_ref(self):
        from mod.reference import RefAnalysis
        return RefAnalysis(self)

    def init_eff(self):
        from mod.efficiency import Efficiency
        return Efficiency(self)

    def get_entries(self):
        return self.Data['Tracks']['NTracks'].size

    def load_data(self, test=False):
        if not test:
            self.Converter.run()
            try:
                f = self.load_file()
                # check if data is complete
                _ = f['Tracks']
                _ = f[self.Plane.get_name()]
                return f
            except (KeyError, OSError):
                warning('could not load data file {} -> start with dummy'.format(self.Run.FileName))
        return self.Dummy.load_file()

    def load_file(self):
        return h5py.File(self.Run.FileName, 'r')

    def reload_data(self):
        self.Data = self.load_data()

    def show_structure(self):
        print('.')
        for i0, (key, grp) in enumerate(self.Data.items()):
            if 'Plane' not in key or key == 'Plane0' or key == f'Plane{self.Plane.Number}':
                print(f'├── {key}')
                for k, g in grp.items():
                    print(f'{"│" if i0 < len(self.Data.keys()) - 1 else " "}   ├── {k}')
                    if hasattr(g, 'keys'):
                        for i in g.keys():
                            print(f'│   │   ├── {i}')
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    # todo: add extra class
    def init_cuts(self, redo=False):
        self.Cuts.set_config(self.BeamTest.Tag, self.DUT.Name)
        self.Cuts.register('fid', self.make_fiducial(_redo=redo), 10, 'fid cut')
        self.Cuts.register('mask', self.make_mask(), 20, 'mask pixels')
        self.Cuts.register('charge', self.get_charges(cut=False) != 0, 30, 'events with non-zero charge')
        self.Cuts.register('cluster', self.get_data('Clusters', 'Size', cut=False) > 0, 90, 'tracks with a cluster')

    def add_cuts(self, redo=False):
        self.Cuts.register('res', self.REF.make_dut_residuals(redo), 69, 'small residuals to REF plane')
        self.Cuts.register('triggerphase', self.make_trigger_phase(redo=redo), 61, 'trigger phase')
        self.Cuts.register('tstart', self.make_start_time(redo=redo), 40, 'exclude first events')
        self.Cuts.register('chi2', self.make_chi2(_redo=redo), 50, 'small chi2')

    def add_track_cuts(self, redo=False):
        self.Tracks.Cuts.register('triggerphase', self.make_trigger_phase(tracks=True, redo=redo), 10, 'track trigger phase')
        self.Tracks.Cuts.register('res', self.REF.make_residuals(redo=redo), 20, 'tracks with a small residual in the REF')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, _redo=redo), 30, 'tracks in fiducial area')
        self.Tracks.Cuts.register('tstart', self.make_start_time(tracks=True, redo=redo), 40, 'exclude first events')
        self.Tracks.Cuts.register('chi2', self.make_chi2(tracks=True, _redo=redo), 50, 'small chi2')

    def reload_cuts(self, redo=False):
        self.init_cuts(redo)
        self.add_cuts(redo)
        self.add_track_cuts(redo)

    def activate_surface(self, on=True):
        self.Cuts.register('fid', self.make_fiducial(surface=on), 10, 'fid cut')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, surface=on), 30)
        self.Surface = on

    def deactivate_surface(self):
        self.activate_surface(False)

    def make_start_time(self, tracks=False, redo=False):
        def f():
            t0 = self.Cuts.get_config('start time', dtype=int) * 60
            t = self.get_time(cut=False, trk_cut=-1 if not tracks else False)
            return t > t[0] + t0
        return array(do_hdf5(self.make_hdf5_path('time', sub_dir='tracks' if tracks else None), f, redo))

    def make_res(self):
        x, y = self.get_du(cut=0), self.get_dv(cut=0)
        x0, y0 = self.get_du(), self.get_dv()
        mx, my = mean(x0[abs(x0) < .3]), mean(y0[abs(y0) < .3])
        return sqrt((x-mx) ** 2 + (y-my) ** 2) < .1

    @staticmethod
    def p_in_poly(p, poly: TCutG):
        return poly.IsInside(*p)

    @save_hdf5('Fid', sub_dir='cuts', arr=True, dtype='?', suf_args='all')
    @parallel('p_in_poly', 'fiducial cut')
    def make_fiducial(self, tracks=False, surface=False, _redo=False):
        x, y = self.Tracks.get_coods(local=True, trk_cut=False) if tracks else self.get_xy(local=True, cut=False)
        return array([x, y]).T, self.get_fid(off=0, surface=surface)

    @save_hdf5('Chi2', sub_dir='cuts', arr=True, suf_args='all')
    def make_chi2(self, tracks=False, q=None, _redo=False):
        x = self.Tracks.get_chi2(cut=False, trk_cut=False if tracks else -1)
        return x < quantile(x, choose(q, self.Cuts.get_config('chi2 quantile', dtype=float)))

    def make_mask(self):
        x, y = self.get_xy(local=True, cut=False)
        mx, my = array(self.Cuts.get_config('mask')).T
        return all([invert((x >= mx[i] - .5) & (x <= mx[i] + .5) & (y >= my[i] - .5) & (y <= my[i] + .5)) for i in range(mx.size)], axis=0)

    def get_fid(self, off=-.5, surface=None, **dkw):
        x, y = self.Cuts.get_fid_config(choose(surface, self.Surface)) + off
        return self.Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=f'fid{self.Surface:d}'))

    def get_full_size(self, off=0, **dkw):
        x, y = make_box_args(*self.Cuts.get_config('full size')[[0, 2, 1, 3]]) + off
        return self.Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=f'full size'))

    def make_trigger_phase(self, tracks=False, redo=False):
        def f():
            tp = self.get_trigger_phase(cut=False, trk_cut=False if tracks else -1)
            low, high = self.Cuts.get_config('trigger phase')
            return (tp >= low) & (tp <= high)
        return array(do_hdf5(self.make_hdf5_path('tp', sub_dir='tracks' if tracks else None), f, redo))

    def remove_metadata(self):
        for p in self.get_meta_files():
            remove_file(p)
        self.Cuts.set_config(self.BeamTest.Tag, self.DUT.Name)
        self.init_cuts()
        self.add_cuts()
        self.add_track_cuts()

    def ev2dut(self, x):
        return self.ev2trk(x)[self.Cuts['cluster']]

    def dut2ev(self, x):
        c0 = self.Cuts['cluster'].copy()
        c0[c0] = x
        c = zeros(self.NEvents, '?')
        c[self.Tracks.get_events()[c0]] = True
        return c

    def ev2trk(self, x):
        return x[self.Tracks.get_events()]

    def trk2ev(self, x):
        c = zeros(self.NEvents, '?')
        c[self.Tracks.get_events()[x]] = True
        return c
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region DATA
    def get_time(self, cut=None, trk_cut: Any = -1):
        t = array(self.Data['Event']['Time']).astype('f8') + self.Run.StartTime
        t = t[self.Tracks.get_events()]  # map the event times to the track times
        return t[self.Cuts.get('cluster')()][self.Cuts(cut)] if trk_cut == -1 else t[self.Tracks.Cuts(trk_cut)]

    def get(self):
        return self.Data[f'Plane{self.Plane.Number}']

    def get_group(self, grp):
        return self.get()[grp]

    def get_data(self, grp, key=None, cut=None):
        data = self.get_group(grp)
        data = array(data) if key is None else array(data[key])
        return data[self.Cuts(cut)]

    def get_charges(self, e=False, cut=None):
        return self.get_data('Clusters', 'Charge', cut) * (self.DUT.VcalToEl if e else 1)

    def get_charge(self, cut=None):
        values = self.get_charges(cut=cut)
        m, s = mean_sigma(values)
        return ufloat(m, s / sqrt(values.size))

    def get_x(self, cut=None):
        return self.get_data('Clusters', 'X', cut)

    def get_y(self, cut=None):
        return self.get_data('Clusters', 'Y', cut)

    def get_u(self, cut=None):
        return self.get_data('Clusters', 'U', cut=cut)

    def get_v(self, cut=None):
        return self.get_data('Clusters', 'V', cut=cut)

    def get_du(self, cut=None, x=None):
        return self.get_u(cut) - choose(x, self.Tracks.get_u, cut=cut)

    def get_dx(self, cut=None, x=None):
        return self.get_x(cut) - choose(x, self.Tracks.get_x, cut=cut)

    def get_dv(self, cut=None, x=None):
        return self.get_v(cut=cut) - choose(x, self.Tracks.get_v, cut=cut)

    def get_dy(self, cut=None, x=None):
        return self.get_y(cut) - choose(x, self.Tracks.get_y, cut=cut)

    def get_xy(self, local=True, cut=None):
        return array([self.get_x(cut=cut) if local else self.get_u(cut), self.get_y(cut=cut) if local else self.get_v(cut)])

    def get_residuals(self, cut=None):
        return sqrt(self.get_du(cut) ** 2 + self.get_dv(cut) ** 2)
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_start_time(self):
        return datetime.fromtimestamp(self.Run.StartTime)

    def get_end_time(self):
        return datetime.fromtimestamp(self.Run.EndTime)

    def get_mask(self):
        return self.get_data('Mask', cut=False)

    def t_args(self, rel_t=False):
        return {'x_tit': 'Time [hh:mm]', 't_ax_off': self.Run.StartTime if rel_t else 0}

    @staticmethod
    def get_ax_tits(local=True):
        return {'x_tit': 'Column', 'y_tit': 'Row'} if local else {'x_tit': 'X [mm]', 'y_tit': 'Y [mm]'}

    def get_track_data(self, grp, key=None, cut=None, trk_cut=-1):
        data = self.get_data(grp, key, cut=False)
        return data[self.Tracks.Cuts(trk_cut)] if type(trk_cut) == ndarray or trk_cut != -1 else data[self.Cuts.get('cluster').Values][self.Cuts(cut)]

    def get_cluster_size(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('Clusters', 'Size', cut=cut, trk_cut=trk_cut)

    def get_trigger_phase(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('TriggerPhase', cut=cut, trk_cut=trk_cut)

    def get_ev(self, cut=None):
        """ returns: event numbers of the dut cut. """
        return where(self.dut2ev(self.Cuts(cut)))[0]

    def get_events(self):
        """ returns: event numbers with clusters. """
        return self.Tracks.get_events(self.Cuts.get('cluster'))

    def get_segments(self, nx, ny, width=False):
        x0, x1, y0, y1 = self.Cuts.get_config('full size')
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

    @property
    def loc_bins(self):
        return bins.get_local(self.Plane)

    @property
    def ph_tit(self):
        return 'Pulse Height [vcal]'

    def get_ph_tit(self, e=False):
        return 'Charge [e]' if e else self.ph_tit
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_mask(self, **dkw):
        self.Draw.histo_2d(*self.get_mask().T, self.loc_bins, 'Masked Pixels', **prep_kw(dkw, **self.get_ax_tits(), fill_color=1, rm=.03, draw_opt='box'))

    def draw_occupancy(self, local=True, bw=1, cut=None, fid=False, **dkw):
        x, y = self.get_xy(local, self.Cuts.get_nofid(cut, fid))
        title = f'{"Local" if local else "Global"} Cluster Occupancy'
        self.Draw.histo_2d(x, y, bins.get_xy(local, self.Plane, bw), title, **prep_kw(dkw, x_tit='Column', y_tit='Row', stats=set_statbox(entries=True, m=True)))

    def draw_hit_map(self, res=.3, local=True, cut=None, trk_cut=-1, fid=False, **dkw):
        x, y = self.Tracks.get_coods(local, self.Cuts.get_nofid(cut, fid), trk_cut)
        self.Draw.histo_2d(x, y, bins.get_xy(local, self.Plane, res), 'Hit Map', **prep_kw(dkw, **self.get_ax_tits(local)))

    def draw_cluster_size(self, cut=None, trk_cut=-1, **dkw):
        v = self.get_cluster_size(cut, trk_cut)
        self.Draw.distribution(v, bins.make(0, 10), 'Cluster Size', **prep_kw(dkw, x_tit='Cluster Size'))

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        cut = self.Cuts.get_nofid(cut, fid)
        x, y = self.Tracks.get_coods(local, cut)
        self.Draw.prof2d(x, y, self.get_cluster_size(cut), bins.get_xy(local, self.Plane, res), 'Cluster Size', **prep_kw(dkw, z_tit='Cluster Size', **self.get_ax_tits(local)))

    def draw_trigger_phase(self, cut=None, trk_cut=-1):
        cut, trk_cut = self.Cuts.exclude('triggerphase', cut), self.Tracks.Cuts.exclude('triggerphase', trk_cut) if trk_cut != -1 else trk_cut
        h = self.Draw.distribution(self.get_trigger_phase(cut, trk_cut), bins.get_triggerphase(), 'Trigger Phase', x_tit='Trigger Phase', y_off=1.8, lm=.13)
        format_histo(h, y_range=[0, h.GetMaximum() * 1.1])
        update_canvas()

    def draw_time(self, show=True):
        t = self.get_time()
        g = self.Draw.profile(arange(t.size), t, bins.make(0, t.size, sqrt(t.size)), x_tit='Event Number', y_tit='Time [hh:mm]', show=show, draw_opt='al', graph=True)
        set_time_axis(g, axis='Y')

    def draw_grid(self, nx=2, ny=3, w=1, width=False):
        self.Draw.grid(*self.get_segments(nx, ny, width), w)

    def draw_x(self, res=.1, cluster=False, cut=None):
        x = self.get_x(cut) if cluster else self.Tracks.get_x(cut)
        self.Draw.distribution(x, bins.make(0, self.Plane.NCols, res, last=True))

    def draw_inpixel_map(self, res=.1, cut=None, cell=False, show=True):
        x, y = self.expand_inpixel(cell=cell, *self.Tracks.get_coods(cut=cut))
        self.Draw.histo_2d(x, y, bins.get_pixel(self.Plane, res, cell=cell), 'Hit Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region RESIDUALS
    def draw_x_residuals(self, cut=None, **dkw):
        self.Draw.distribution(self.get_du(cut) * 1e3, **prep_kw(dkw, r=[-300, 300], title='X Residuals', x_tit='Residual [#mum]'))

    def draw_y_residuals(self, cut=None, **dkw):
        self.Draw.distribution(self.get_dv(cut) * 1e3, **prep_kw(dkw, r=[-300, 300], title='Y Residuals', x_tit='Residual [#mum]'))

    def draw_xy_residuals(self, cut=None, **dkw):
        x, y = array([f(self.Cuts.exclude('res', cut)) for f in [self.get_du, self.get_dv]]) * 1e3
        self.Draw.histo_2d(x, y, **prep_kw(dkw, q=.1, title='XY Residual', x_tit='dX [#mum]', y_tit='dY [#mum]'))

    def draw_residuals(self, **dkw):
        self.Draw.distribution(self.get_residuals() * 1e3, **prep_kw(dkw, title='Residuals', x_tit='Residual [#mum]'))

    def draw_residuals_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        (x, y), z_ = [f(cut=self.Cuts.get_nofid(cut, fid)) for f in [partial(self.Tracks.get_coods, local=local), self.get_residuals]]
        self.Draw.prof2d(x, y, z_ * 1e3, bins.get_xy(local, self.Plane, res), 'Residuals', **prep_kw(dkw, z_tit='Residuals [#mum]', **self.get_ax_tits(local)))

    def draw_angle(self, x, y, prof=False, xb=True, local=False, **dkw):
        b = (bins.get_x if xb else bins.get_y)(self.Plane, local=local) + find_bins(y)
        return arctan(FitRes((self.Draw.profile if prof else self.Draw.histo_2d)(x, y, b[:2 if prof else 4], graph=True, **dkw).Fit('pol1', 'qs'))[1].n)

    def draw_udv(self, cut=None, tv=None, prof=False, **dkw):
        return self.draw_angle(self.get_u(cut), self.get_dv(cut, tv), prof, **prep_kw(dkw, title='X dY', x_tit='X [mm]', y_tit='dY [mm]'))

    def draw_vdu(self, cut=None, tu=None, prof=False, **dkw):
        return self.draw_angle(self.get_v(cut), self.get_du(cut, tu), prof, xb=False, **prep_kw(dkw, title='Y dX', x_tit='Y [mm]', y_tit='dX [mm]'))

    def draw_xdy(self, cut=None, ty=None, prof=False, **dkw):
        return self.draw_angle(self.get_x(cut), self.get_dy(cut, ty), prof, local=True, **prep_kw(dkw, title='X dY', x_tit='X [Column]', y_tit='dY [Row]'))

    def draw_ydx(self, cut=None, tx=None, prof=False, **dkw):
        return self.draw_angle(self.get_y(cut), self.get_dx(cut, tx), prof, xb=False, local=True, **prep_kw(dkw, title='Y dX', x_tit='Y [row]', y_tit='dX [col]'))
    # endregion RESIDUALS
    # ----------------------------------------

    # ----------------------------------------
    # region CORRELATION
    def corr_cuts(self, pl=2):
        cdut = self.ev2dut(self.Tel.get_n_clusters(pl, cut=False) == 1)
        return cdut, self.dut2ev(cdut)

    def draw_x_correlation(self, res=1, pl=2, **dkw):
        cdut, ctel = self.corr_cuts(pl)
        x, y = self.get_x(cdut), self.Tel.get_x(pl, cut=ctel)
        self.Draw.histo_2d(x, y, find_bins(x, w=res) + find_bins(y, q=.2), **prep_kw(dkw, title='X Corr', x_tit='X DUT [Column]', y_tit=f'X Plane {pl} [Column]'))

    def draw_y_correlation(self, res=1, pl=2, **dkw):
        cdut, ctel = self.corr_cuts(pl)
        x, y = self.get_y(cdut), self.Tel.get_y(pl, cut=ctel)
        self.Draw.histo_2d(x, y, find_bins(x, w=res) + find_bins(y, q=.2), **prep_kw(dkw, title='Y Corr', x_tit='Y DUT [Column]', y_tit=f'Y Plane {pl} [Column]'))

    def draw_correlation_trend(self, pl=2, thresh=.2, **dkw):
        cdut, ctel = self.corr_cuts(pl)
        d0, d1, t = self.get_xy(cut=cdut), self.Tel.get_xy(pl, cut=ctel), self.get_time(cdut)
        g = [self.Draw.graph(*get_3d_correlations(self.Draw.histo_3d(t, d0[i], d1[i]), thresh=thresh), **self.t_args(), y_tit='Correlation Factor', show=False) for i in range(2)]
        return self.Draw.multigraph(g, 'CorrFac', ['x', 'y'], draw_opt='pl', **prep_kw(dkw, **self.t_args()))

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
    def draw_ped_map(self, c_max, c_min=None):
        charge = self.get_charges(cut=False)
        cut = Cut('charge', charge < c_max if c_min is None else (charge > c_min) & (charge < c_max)) + self.Cuts.exclude('fid')
        self.draw_hit_map(cut=cut)

    def draw_signal_map(self, res=.3, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cuts.get_nofid(cut, fid)) for f in [self.Tracks.get_coods, self.get_charges]]
        self.Draw.prof2d(x, y, z_, bins.get_local(self.Plane, res), 'Charge Map', **prep_kw(dkw, leg=self.get_fid(), z_tit=self.ph_tit, **self.get_ax_tits()))

    def draw_signal_occupancy(self, fid=False, cut=None, **dkw):
        (x, y), z_ = [f(cut=self.Cuts.get_nofid(cut, fid)) for f in [self.get_xy, self.get_charges]]
        self.Draw.prof2d(x, y, z_, self.loc_bins, 'Charge Occupancy', **prep_kw(dkw, leg=self.get_fid(), z_tit=self.ph_tit, **self.get_ax_tits()))

    def draw_signal_distribution(self, cut=None, draw_thresh=False, e=False, **dkw):
        trim = self.Calibration.Trim * (self.DUT.VcalToEl if e else 1)
        legs = [self.Draw.vertical_line(trim, 0, 1e5, w=2, show=False), self.Draw.tlatex(.2, .5, f'Threshold = {trim:.0f} vcal', ndc=True, angle=90, size=.04, show=False)] if draw_thresh else None
        return self.Draw.distribution(self.get_charges(e, cut), **prep_kw(dkw, title='PH', x_tit=self.ph_tit, leg=legs))

    def draw_charge_distribution(self, cut=None, draw_thresh=False, **dkw):
        return self.draw_signal_distribution(cut, draw_thresh, e=True, **dkw)

    def draw_signal_vs_tp(self, cut=None, **dkw):
        x, y = [f(cut=self.Cuts.exclude('triggerphase', cut)) for f in [self.get_trigger_phase, self.get_charges]]
        self.Draw.profile(x, y, bins.get_triggerphase(), 'Charge vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', y_tit=self.ph_tit, graph=True))

    def draw_signal_trend(self, bw=None, e=False, cut=None, **dkw):
        x, y = self.get_time(cut), self.get_charges(e=e, cut=cut)
        g = self.Draw.profile(x, y, find_bins(x, w=bw), **self.t_args(), graph=True, y_tit=self.get_ph_tit(e), stats=True, show=False)
        return self.Draw(g, **prep_kw(dkw, y_range=ax_range(get_graph_y(g, err=False), fl=1, fh=2)))

    def fit_signal(self, bw=None, e=False, **dkw):
        g = self.draw_signal_trend(bw, e, show=False)
        fit = FitRes(g.Fit('pol0', 'sq'))
        self.Draw(g, **prep_kw(dkw, stats=set_statbox(fit=True), show=False))
        return fit

    def draw_inpixel_charge(self, res=.1, cut=None, show=True, cell=False):
        (x, y), c = self.Tracks.get_coods(cut), self.get_charges(cut=cut)
        x, y, c = self.expand_inpixel(x, y, c, cell)
        self.Draw.prof2d(x, y, c, bins.get_pixel(self.Plane, res, cell=cell), 'Charge Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
    # endregion SIGNAL
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

    def estimate_mean(self, n=10, fit_range=None):
        fit_range = [8000, 30000] if fit_range is None else fit_range
        self.fit_langau(show=False, fit_range=fit_range)
        values = []
        self.PBar.start(n)
        for i in range(n):
            h = TH1F('hph{}'.format(i), 'Pulse Height Distribution', 200, -5000, 70000)
            h.FillRandom('Fitfcn_hph', 25000)
            values.append(h.GetMean())
            self.PBar.update(i)
        return mean_sigma(values)
