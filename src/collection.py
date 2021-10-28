from src.dut_analysis import DUTAnalysis
from plotting.draw import *

r = []
for i in range(9):
    print(i)
    r.append(DUTAnalysis('data/Clustered_0{}.root'.format(i)))

draw = Draw()
x = [t.DUT.Bias * ufloat(1, .01) for t in r]
y = [t.get_efficiency() for t in r]

draw.graph(x, y, x_tit='Voltage [V]', y_tit='Charge [vcal]')
