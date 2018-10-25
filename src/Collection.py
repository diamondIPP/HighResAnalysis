from RunAnalysis import RunAnalysis
from draw import *
from uncertainties import ufloat

r = []
for i in xrange(9):
    print i
    r.append(RunAnalysis('data/Clustered_0{}.root'.format(i)))


d = Draw()
y=[]
x=[-60,-70,-60,-50,-40,-30,-20,-10, 0]
x = [make_ufloat([i, -.01 *i]) for i in x]
for t in r:
    print t.RunNumber
    fit = t.fit_efficiency(3*60, cut='Timing<5')
    y.append(make_ufloat([fit.Parameter(0), fit.ParError(0)]))

g = d.make_tgrapherrors('g', 'g', x=x, y=y)
d.format_histo(g, x_tit='Voltage [V]', y_tit='Charge [vcal]', y_off=1.3)
d.draw_histo(g, lm=1.1)