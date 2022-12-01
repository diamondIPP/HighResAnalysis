# --------------------------------------------------------
#       define a DUT
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from utility.utils import critical, ufloat, array, Dir
from plotting.draw import Draw, arange, prep_kw, add_perr, Config
from src.analysis import Analysis


class Device:
    """ parent class with information about a single device. """
    def __init__(self, number=1, name='Name', typ=None, has_ref=False):
        self.Number = number
        self.Name = name
        self.Type = typ
        self.Plane = self.init_plane(has_ref)

    def __str__(self):
        return self.Name

    def __repr__(self):
        return f'{self.Type} {self.Number}, {self}'

    def init_plane(self, has_ref):
        return Plane(Analysis.Config.getint('TELESCOPE', 'planes') + self.Number + int(has_ref), self.Type)


class REF(Device):
    """ Class with information about the reference plane. """
    def __init__(self, number=0, name='REF'):
        super().__init__(number, name, typ='REF')


class DUT(Device):
    """ Class with all information about a single DUT. """
    def __init__(self, number=1, run_log: dict = None, has_ref=False):

        # Info
        super().__init__(number, run_log['duts'][number], typ='DUT', has_ref=has_ref)
        self.Bias = int(run_log[f'hv'][self.Number])
        self.Position = int(run_log[f'dut position'][self.Number])

        # Specs
        self.Info = self.load_specs()
        self.Irradiation = self.Info.get_value('irradiation', default={})
        self.Thickness = self.Info.get_value('thickness', default=500)
        self.CCD = self.Info.get_value('CCD')
        self.Size = self.Info.get_value('size', default=[5, 5])
        self.Cells = self.Info.get_value('cells')
        if self.Cells is not None:
            self.NColumns = 2 * self.Cells[0] * self.Cells[1] + sum(self.Cells) + 1
            self.ColumnDiameter = add_perr(self.Info.get_float('column diameter'), .05)
            self.PXY = array(self.Info.get_list('cell size'))
            self.PXYu = self.PXY * 1e3  # in um
            self.PX, self.PY = self.PXY
            self.PXu, self.PYu = self.PXYu
        self.VcalToEl = Analysis.Config.get_float('DUT', 'vcal to electrons')

    def __repr__(self):
        return f'{super().__repr__()}, Bias: {self.Bias:1.0f}V'

    def load_specs(self):
        f = Dir.joinpath('config', 'dia_info.json')
        return Config(f, section=self.Name, from_json=True)

    def get_irradiation(self, tc):
        return self.Irradiation[tc] if tc in self.Irradiation else critical('Please add "{}" to the irradiation file for {}'.format(self.Name, tc))

    def load_spec(self, section, typ=None, lst=False, error=None, default=None):
        spec = default if section not in self.Info or self.Info[section] == 'None' else self.Info[section] if typ is None else typ(self.Info[section])
        return [v if error is None else ufloat(v, error) for v in spec] if lst and spec is not None else ufloat(spec, error) if error is not None and spec is not None else spec

    def set_number(self, value):
        self.Number = value


class Plane:
    """ Class with all information about a single pixel plane. """
    def __init__(self, n, typ='DUT', rotated=False):

        config = Analysis.Config(typ)
        self.IsDUT = 'DUT' in config.Section
        self.Number = n
        self.Type = config.get_value('name')
        self.NCols, self.NRows = config.get_value('pixel')
        self.NPixels = self.NCols * self.NRows
        self.PXY = array(config.get_value('pitch'))
        self.PX, self.PY = self.PXY
        self.R = self.PX / self.PY
        self.M = array([[self.PX, 0], [0, self.PY]])
        self.W, self.H = self.PX * self.NCols, self.PY * self.NRows
        self.Rotated = rotated

    def __str__(self):
        return f'Plane{self.Number}'

    def __repr__(self):
        return f'{self}, {self.Type.upper()}: {self.NCols}x{self.NRows} pixels ({self.PX * 1e3:1.1f}x{self.PY * 1e3:1.1f}μm)'

    def __add__(self, other):
        self.Number += other
        return self

    def get_max_width(self):
        return max(self.get_x_width(), self.get_y_width())

    def get_x_width(self):
        return self.PX * self.NCols

    def get_y_width(self):
        return self.PY * self.NRows

    def get_grid(self, off=-.5, **dkw):
        return Draw.grid(arange(self.NCols + 1) + off, arange(self.NRows + 1) + off, **prep_kw(dkw, show=False))

    def draw_grid(self, off=-.5, **dkw):
        self.get_grid(off, **prep_kw(dkw, show=True))
