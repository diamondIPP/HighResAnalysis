#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a single DUT
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F
from numpy import cumsum, split, histogram, invert, histogram2d, max, arctan

from src.analysis import *
from src.calibration import Calibration
from src.cern_run import CERNRun
from src.converter import Converter
from src.currents import Currents
from src.cut import Cut
from src.desy_converter import DESYConverter
from src.desy_run import DESYRun
from src.dummy import Dummy
from plotting.fit import *
from src.reference import RefAnalysis
from src.utils import *
from src.telescope import TelescopeAnalysis
from src.tracks import TrackAnalysis
import src.bins as bins


# todo add ADC to data!
# TODO: add residuals to data!
class DUTAnalysis(Analysis):

    def __init__(self, run_number, dut_number, test_campaign, single_mode=True, verbose=True, test=False):

        Analysis.__init__(self, test_campaign, meta_sub_dir='DUT', verbose=verbose)
        self.print_start(run_number)

        # MAIN
        self.Run = self.run(run_number, dut_number, self.BeamTest.Path, self.Config, single_mode)
        self.DUT = self.Run.DUT
        self.Plane = self.DUT.Plane

        # DATA
        self.Converter = self.converter(self.BeamTest.Path, self.Run.Number, self.Config)
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
        self.Telescope = TelescopeAnalysis(self)
        self.Tracks = TrackAnalysis(self)
        self.REF = RefAnalysis(self)
        self.Currents = Currents(self)
        self.add_cuts()
        self.add_track_cuts()

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
        self.Cuts.register('fid', self.make_fiducial(redo=redo), 10, 'fid cut')
        self.Cuts.register('mask', self.make_mask(), 20, 'mask pixels')
        self.Cuts.register('charge', self.get_charges(cut=False) != 0, 30, 'events with non-zero charge')
        self.Cuts.register('cluster', self.get_data('Clusters', 'Size', cut=False) > 0, 90, 'tracks with a cluster')

    def add_cuts(self, redo=False):
        self.Cuts.register('res', self.REF.make_dut_residuals(redo), 69, 'small residuals to REF plane')
        self.Cuts.register('triggerphase', self.make_trigger_phase(redo=redo), 61, 'trigger phase')
        self.Cuts.register('tstart', self.make_start_time(redo=redo), 40, 'exclude first events')
        self.Cuts.register('chi2', self.make_chi2(redo=redo), 50, 'small chi2')

    def add_track_cuts(self, redo=False):
        self.Tracks.Cuts.register('triggerphase', self.make_trigger_phase(tracks=True, redo=redo), 10, 'track trigger phase')
        self.Tracks.Cuts.register('res', self.REF.make_residuals(redo=redo), 20, 'tracks with a small residual in the REF')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, redo=redo), 30, 'tracks in fiducial area')
        self.Tracks.Cuts.register('tstart', self.make_start_time(tracks=True, redo=redo), 40, 'exclude first events')
        self.Tracks.Cuts.register('chi2', self.make_chi2(tracks=True, redo=redo), 50, 'small chi2')

    def reload_cuts(self, redo=False):
        self.init_cuts(redo)
        self.add_cuts(redo)
        self.add_track_cuts(redo)

    def activate_surface(self, on=True):
        option = '{}fiducial'.format('surface ' if on else '')
        self.Cuts.register('fid', self.make_fiducial(option=option), 10, 'fid cut')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, option=option), 30)
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
        x, y = self.Tracks.get_coods(local=True, trk_cut=False) if tracks else self.get_coods(local=True, cut=False)
        return array([x, y]).T, self.get_fid(off=0, surface=surface)

    def make_chi2(self, tracks=False, redo=False):
        def f():
            chi2 = self.Tracks.get_chi2(cut=False, trk_cut=-1)
            chi_max = quantile(chi2, self.Cuts.get_config('chi2 quantile', dtype=float))
            return self.Tracks.get_chi2(cut=False, trk_cut=False if tracks else -1) < chi_max
        return array(do_hdf5(self.make_hdf5_path('chi2', sub_dir='tracks' if tracks else None), f, redo))

    def make_mask(self):
        x, y = self.get_coods(local=True, cut=False)
        mx, my = array(self.Cuts.get_config('mask')).T
        return all([invert((x >= mx[i] - .5) & (x <= mx[i] + .5) & (y >= my[i] - .5) & (y <= my[i] + .5)) for i in range(mx.size)], axis=0)

    def get_fid(self, off=-.5, surface=None, **dkw):
        p = self.Cuts.get_fid_config(choose(surface, self.Surface)) + off
        x, y = make_box_args(*p[[0, 2, 1, 3]] + array([0, 0, 1, 1])) if p.size == 4 else p
        return self.Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=f'fid{self.Surface:d}'))

    def make_trigger_phase(self, tracks=False, redo=False):
        def f():
            tp = self.get_trigger_phase(cut=False, trk_cut=False if tracks else -1)
            low, high = self.Cuts.get_config('trigger phase')
            return (tp >= low) & (tp <= high)
        return array(do_hdf5(self.make_hdf5_path('tp', sub_dir='tracks' if tracks else None), f, redo))

    def make_correlation(self, plane=2):
        n = self.Telescope.get_n('Clusters', plane, cut=False)
        return (n[self.Tracks.get_events()] == 1)[self.Cuts.get('cluster')()]

    def remove_metadata(self):
        for p in self.get_meta_files():
            remove_file(p)
        self.Cuts.set_config(self.BeamTest.Tag, self.DUT.Name)
        self.init_cuts()
        self.add_cuts()
        self.add_track_cuts()
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_start_time(self):
        return datetime.fromtimestamp(self.Run.StartTime)

    def get_end_time(self):
        return datetime.fromtimestamp(self.Run.EndTime)

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

    def get_residuals(self, cut=None):
        return sqrt(self.get_du(cut) ** 2 + self.get_dv(cut) ** 2)

    def get_coods(self, local=True, cut=None):
        return array([self.get_x(cut=cut) if local else self.get_u(cut), self.get_y(cut=cut) if local else self.get_v(cut)])

    def get_mask(self):
        return self.get_data('Mask', cut=False)

    def get_time_args(self, rel_t=False):
        return {'x_tit': 'Time [hh:mm]', 't_ax_off': self.Run.StartTime if rel_t else 0}

    @staticmethod
    def get_ax_tits(local=True):
        return {'x_tit': 'Column', 'y_tit': 'Row'} if local else {'x_tit': 'X [mm]', 'y_tit': 'Y [mm]'}

    def get_track_data(self, grp, key=None, cut=None, trk_cut=-1):
        data = self.get_data(grp, key, cut=False)
        return data[self.Tracks.Cuts(trk_cut)] if type(trk_cut) == ndarray or trk_cut != -1 else data[self.Cuts.get('cluster').Values][self.Cuts(cut)]

    def get_cluster_size(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('Clusters', 'Size', cut=cut, trk_cut=trk_cut)

    def get_efficiencies(self, trk_cut=None):
        return (self.get_cluster_size(trk_cut=trk_cut) > 0).astype('u2') * 100

    def get_segment_efficiecies(self, nx=2, ny=3, cut=None):
        return get_2d_hist_vec(self.draw_efficiency_map(local=True, cut=cut, binning=bins.make2d(*self.get_segments(nx, ny)), show=False), err=False, flat=False)

    def get_efficiency(self, trk_cut=None, prnt=True):
        eff = calc_eff(values=self.get_efficiencies(trk_cut))
        self.info('{:.1f}(-{:.1f}+{:.1f})%'.format(*eff), prnt=prnt)
        return eff

    def get_trigger_phase(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('TriggerPhase', cut=cut, trk_cut=trk_cut)

    def get_events(self):
        """ returns: event numbers with clusters. """
        return self.Tracks.get_events(self.Cuts.get('cluster'))

    def get_segments(self, nx, ny):
        x0, x1, y0, y1 = self.Cuts.get_config('full size', lst=True)
        return arange(x0, x1 + bool((x1 - x0 + 1) % nx) * nx + 1.1, nx, dtype='u2'), arange(y0, y1 + bool((y1 - y0 + 1) % ny) * ny + 1.1, ny, dtype='u2')

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
        x, y = self.get_coods(local, self.Cuts.get_nofid(cut, fid))
        title = f'{"Local" if local else "Global"} Cluster Occupancy'
        self.Draw.histo_2d(x, y, bins.get_coods(local, self.Plane, bw), title, **prep_kw(dkw, x_tit='Column', y_tit='Row', stats=set_statbox(entries=True, m=True)))

    def draw_hit_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        x, y = self.Tracks.get_coods(local, self.Cuts.get_nofid(cut, fid))
        self.Draw.histo_2d(x, y, bins.get_coods(local, self.Plane, res), 'Hit Map', **prep_kw(dkw, **self.get_ax_tits(local)))

    def draw_cluster_size(self, cut=None, trk_cut=-1, **dkw):
        v = self.get_cluster_size(cut, trk_cut)
        self.Draw.distribution(v, bins.make(0, 10), 'Cluster Size', **prep_kw(dkw, x_tit='Cluster Size'))

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, **dkw):
        cut = self.Cuts.get_nofid(cut, fid)
        x, y = self.Tracks.get_coods(local, cut)
        self.Draw.prof2d(x, y, self.get_cluster_size(cut), bins.get_coods(local, self.Plane, res), 'Cluster Size', **prep_kw(dkw, z_tit='Cluster Size', **self.get_ax_tits(local)))

    def draw_trigger_phase(self, cut=None, trk_cut=-1):
        cut, trk_cut = self.Cuts.exclude('triggerphase', cut), self.Tracks.Cuts.exclude('triggerphase', trk_cut) if trk_cut != -1 else trk_cut
        h = self.Draw.distribution(self.get_trigger_phase(cut, trk_cut), bins.get_triggerphase(), 'Trigger Phase', x_tit='Trigger Phase', y_off=1.8, lm=.13)
        format_histo(h, y_range=[0, h.GetMaximum() * 1.1])
        update_canvas()

    def draw_time(self, show=True):
        t = self.get_time()
        g = self.Draw.profile(arange(t.size), t, bins.make(0, t.size, sqrt(t.size)), x_tit='Event Number', y_tit='Time [hh:mm]', show=show, draw_opt='al', graph=True)
        set_time_axis(g, axis='Y')

    def draw_segments(self, nx=2, ny=3, w=1):
        x, y = array(self.get_segments(nx, ny)) - .5
        self.Draw.grid(x, y, w)

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
        self.Draw.prof2d(x, y, z_ * 1e3, bins.get_coods(local, self.Plane, res), 'Residuals', **prep_kw(dkw, z_tit='Residuals [#mum]', **self.get_ax_tits(local)))

    def draw_angle(self, x, y, prof=False, xb=True, local=False, **dkw):
        b = (bins.x if xb else bins.y)(self.Plane, local=local) + find_bins(y)
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
    def draw_correlation(self, mode='y', res=1, plane=2, thresh=.1, show=True):
        return self.Telescope.draw_dut_correlation(mode, res, plane, thresh, show)

    def draw_x_correlation(self, res=1, plane=2, show=True):
        self.draw_correlation('x', res, plane, show)

    def draw_y_correlation(self, res=1, plane=2, show=True):
        self.draw_correlation('y', res, plane, show)

    def draw_correlation_trend(self, plane=2, bw=120, thresh=.5, **dkw):
        v0, v1 = self.get_y(cut=self.make_correlation(plane)), self.Telescope.get_y(plane, cut=self.Telescope.make_correlation(plane))
        splits, t = histogram(self.get_time(self.make_correlation(plane)), bins.get_time(self.get_time(), min(int(self.Duration / 3), bw))[1])
        hs = [histogram2d(v[0], v[1], [bins.get_local_y(self.Plane, 2)[1], bins.get_local_y(self.Telescope.Plane, 4)[1]])[0] for v in split(array([v0, v1]), cumsum(splits)[:-1], axis=1)]
        # take only bins with >25% of max entries (noise filter)
        values = [corrcoef(array(where(h > max(h) * thresh)).repeat(h[h > max(h) * thresh].astype('i'), axis=1).T, rowvar=False)[0][1] for h in hs]
        return self.Draw.graph(t[:-1] + diff(t) / 2, y=values, **prep_kw(dkw, y_tit='Correlation Factor', **self.get_time_args(), y_range=[0, 1.05], draw_opt='ap'))

    def draw_alignment(self, plane=2, bw=120, thresh=.5, show=True):
        x, c = get_graph_vecs(self.draw_correlation_trend(plane, bw, show=False), err=False)
        r = [1 if ic < thresh else 2 for ic in c]
        x, y = x.repeat(r), ones(sum(r))
        binning = bins.get_time(self.get_time(), min(int(self.Duration / 3), bw)) + [3, 0, 3]
        gStyle.SetPalette(3, array([1, 2, 3], 'i'))
        self.Draw.histo_2d(x, y, binning, 'Event Alignment', show=show, t_ax_off=0, x_tit='Time [hh:mm]', y_tit='Alignment', stats=False, l_off_y=99, center_y=True, draw_opt='col', rm=.03)
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
        (x, y), z_ = [f(cut=self.Cuts.get_nofid(cut, fid)) for f in [self.get_coods, self.get_charges]]
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
        g = self.Draw.profile(x, y, find_bins(x, w=bw), **self.get_time_args(), graph=True, y_tit=self.get_ph_tit(e), stats=True, show=False)
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

    # ----------------------------------------
    # region EFFICIENCY
    def draw_efficiency(self, bin_width=30, show=True):
        t, e = self.get_time(trk_cut=None), self.get_efficiencies()
        p = self.Draw.efficiency(t, e, bins.get_time(t, bin_width), 'Efficiency', x_tit='Time [hh:mm]', y_tit='Efficiency [%]', t_ax_off=0, y_range=[0, 105], show=show, stats=0)
        return p

    def draw_efficiency_vs_trigger_phase(self, show=True):
        cut = self.Tracks.Cuts.exclude('triggerphase')
        x, y = self.get_trigger_phase(trk_cut=cut), self.get_efficiencies(cut)
        return self.Draw.profile(x, y, bins.get_triggerphase(), 'Efficiency vs. Trigger Phase', x_tit='Trigger Phase', y_tit='Efficiency [%]', y_range=[0, 105], show=show)

    def draw_efficiency_map(self, res=.25, local=True, eff=True, both=False, fid=False, cut=None, binning=None, show=True):
        mcut = self.Tracks.Cuts(cut) if fid else self.Tracks.Cuts.exclude('fid')
        x, y = self.Tracks.get_x(trk_cut=mcut, local=local), self.Tracks.get_y(trk_cut=mcut, local=local)
        binning = choose(binning, bins.get_coods, 'None', local, self.Plane, res)
        p = self.Draw.prof2d(x, y, self.get_efficiencies(mcut), binning, 'Efficiency Map', show=show, draw_opt='colz', **self.get_ax_tits(local))
        self.draw_fid(show=not fid and show)
        self.draw_eff_text(self.Surface, cut, eff)
        self.draw_eff_text(not self.Surface, cut, eff and both)
        return p

    def draw_eff_text(self, surface, cut, show=True):
        if show:
            self.activate_surface(surface)
            self.draw_fid()
            x0, x1, y0, y1 = self.Cuts.get_fid_config(surface)
            self.Draw.tlatex(x0 + (x1 - x0) / 2, y0 + (y1 - y0) / 2, '{:2.1f}%'.format(self.get_efficiency(cut, False)[0]), 'eff', 22, size=.04)

    def draw_segment_efficiencies(self, res=.25, local=True, nx=2, ny=3, cut=None, show=True):
        e = self.get_segment_efficiecies(nx, ny, cut)
        self.draw_efficiency_map(res, local, cut=cut, show=show, eff=False)
        x, y = self.get_segments(nx, ny)
        for i, ix in enumerate(x[:-1]):
            for j, iy in enumerate(y[:-1]):
                self.Draw.tlatex(ix + nx / 2 - .5, iy + ny / 2 - .5, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02)
        self.draw_segments(nx, ny)

    def draw_segment_distribution(self, nx=2, ny=3, zsup=False, cut=None, bin_width=None, segments=True, show=True):
        e = self.get_segment_efficiecies(nx, ny, cut).flatten() if segments else get_2d_hist_vec(self.draw_efficiency_map(.5, show=False), err=False)
        binning = bins.make(0, 101.5, choose(bin_width, 1)) if zsup else bins.make(95, 100.5, choose(bin_width, 5.5 / sqrt(e.size) / 2), last=True)
        self.Draw.distribution(e, binning, 'Segment Efficiencies', x_tit='Efficiency [%]', show=show)

    def draw_inpixel_eff(self, res=.1, cut=None, show=True, cell=False):
        (x, y), e = self.Tracks.get_coods(trk_cut=cut), self.get_efficiencies(cut)
        x, y, e = self.expand_inpixel(x, y, e, cell)
        self.Draw.prof2d(x, y, e, bins.get_pixel(self.Plane, res, cell=cell), 'Efficiency Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
    # endregion EFFICIENCY
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
