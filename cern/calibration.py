# --------------------------------------------------------
#       pulse height calibration addon for CERN data
# created on July 2nd 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.calibration import Calibration, Run, critical
from utility.utils import remove_letters


class CERNCalibration(Calibration):

    def __init__(self, run: Run):
        super().__init__(run)

    def load_raw_filename(self):
        files = sorted(list(self.Dir.glob('phCal[0-9]*.dat')), key=lambda x: int(remove_letters(x.name)))
        f = files[next((i - 1 for i, f in enumerate(files) if self.Run.Number < int(remove_letters(f.name.split('-')[0]))), -1)]
        return f if f.exists() else critical(f'could not find adc calibration file {f} ...')

    def get_trim_number(self, n=None):
        f = self.load_raw_filename().stem.split('-')
        return int(f[1]) if len(f) == 2 else None, n
