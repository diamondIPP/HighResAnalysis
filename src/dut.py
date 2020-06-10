# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from analysis import load_json, get_base_dir, OrderedDict, critical, join, ufloat, expanduser, choose
from json import load, loads


class DUT:
    """ Class with all information about a single DUT. """
    def __init__(self, number=1, run_log=None, config=None):

        self.Config = config
        self.Dir = get_base_dir()

        # Info
        self.Number = number
        self.Name = run_log['dut{}'.format(self.Number)]
        self.Bias = run_log['hv{}'.format(self.Number)]
        self.Plane = Plane(config, section='DUT')

        # Specs
        self.Specs = load_json(join(expanduser(self.Config.get('MAIN', 'data directory')), 'dia_info.json'))[self.Name.upper()]
        self.Irradiation = self.load_spec('irradiation')
        self.Thickness = self.load_spec('thickness', typ=int, default=500)
        self.CCD = self.load_spec('CCD', typ=int)
        self.Size = self.load_spec('size', lst=True)
        self.ActiveSize = self.load_spec('active size', lst=True, error=.02)
        self.ActiveArea = self.ActiveSize[0] * self.ActiveSize[1] if self.ActiveSize is not None else None
        self.Pixel = self.load_spec('pixel', lst=True)
        self.NColumns = choose(2 * self.Pixel[0] * self.Pixel[1] + self.Pixel[0] + self.Pixel[1] + 1, default=None, decider=self.Pixel)
        self.ColumnDiameter = self.load_spec('column diameter', typ=float, error=.05)
        self.CellSize = self.load_spec('cell size', typ=int)

    def __str__(self):
        return 'DUT {}, {}, Bias: {:1.0f}V'.format(self.Number, self.Name, self.Bias)

    def __repr__(self):
        return self.__str__()

    def load_irradiation(self):
        with open(join(self.Dir, self.Config.get('MISC', 'irradiation file'))) as f:
            data = load(f)
            return OrderedDict([(key, dic[self.Name]) for key, dic in sorted(data.iteritems()) if self.Name in dic])

    def get_irradiation(self, tc):
        return self.Irradiation[tc] if tc in self.Irradiation else critical('Please add "{}" to the irradiation file for {}'.format(self.Name, tc))

    def load_spec(self, section, typ=None, lst=False, error=None, default=None):
        spec = default if section not in self.Specs or self.Specs[section] == 'None' else self.Specs[section] if typ is None else typ(self.Specs[section])
        return [v if error is None else ufloat(v, error) for v in loads(spec)] if lst and spec is not None else ufloat(spec, error) if error is not None and spec is not None else spec

    def set_number(self, value):
        self.Number = value


class Plane:
    """ Class with all information about a single pixel plane. """
    def __init__(self, config, section='TELESCOPE'):

        self.Name = config.get(section, 'name')
        self.NCols, self.NRows = loads(config.get(section, 'pixel'))
        self.PX, self.PY = loads(config.get(section, 'pitch'))

    def __str__(self):
        return '{} Plane with {}x{} pixels of a size {:1.1f}x{:1.1f}um'.format(self.Name.upper(), self.NCols, self.NRows, self.PX * 1e3, self.PY * 1e3)

    def __repr__(self):
        return self.__str__()
