#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, TLegend, TArrow, TPad, TCutG, TLine, TPaveText, TPaveStats
from ROOT import gROOT, gStyle, kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, kCyan, kTeal
from utils import *
from os.path import dirname, join
from numpy import ndarray, zeros, array
from uncertainties.core import Variable, ufloat
from screeninfo import get_monitors


resolution = None


class Draw:

    def __init__(self, verbose=True, titles=True, save_dir=''):

        # BASICS
        self.Verbose = verbose
        self.Res = load_resolution()
        self.ResultsDir = 'Results'
        self.SaveDir = save_dir
        self.FileTypes = ['root', 'png', 'pdf']

        # COLORS/SETTINGS
        self.count = 0
        self.colors = create_colorlist()
        self.FillColor = 821
        gStyle.SetLegendFont(42)
        self.ActivateTitle = titles
        gStyle.SetOptTitle(titles)

        self.Objects = []

    # ----------------------------------------
    # region BASIC
    def set_save_directory(self, name):
        self.ResultsDir = name

    def make_bias_string(self, bias=None):
        if bias is None:
            return self.make_bias_string(self.bias) if hasattr(self, 'bias') else ''
        pol = 'm' if bias < 0 else 'p'
        return '_{pol}{bias:04d}'.format(pol=pol, bias=int(abs(bias)))

    def get_color(self):
        self.count %= 20
        color = self.colors[self.count]
        self.count += 1
        return color

    def get_colors(self, n):
        return array([self.get_color() for _ in xrange(n)], 'i')

    def reset_colors(self):
        self.count = 0
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    def draw_histo(self, histo, save_name='', show=True, sub_dir=None, lm=.1, rm=.03, bm=.1, tm=None, draw_opt='', x=None, y=None, leg=None, logy=False, logx=False, logz=False,
                   canvas=None, grid=False, gridy=False, gridx=False, prnt=True, phi=None, theta=None):
        return self.save_histo(histo, save_name, show, sub_dir, lm, rm, bm, tm, draw_opt, x, y, leg, logy, logx, logz, canvas, grid, gridx, gridy, False, prnt, phi, theta)

    def draw_axis(self, x1, x2, y1, y2, title, limits=None, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, tick_size=0.03, line=False, opt='+SU', l_off=.01, log=False):
        limits = ([y1, y2] if x1 == x2 else [x1, x2]) if limits is None else limits
        a = TGaxis(x1, y1, x2, y2, limits[0], limits[1], 510, opt + ('G' if log else ''))
        a.SetName(name)
        a.SetLineColor(col)
        a.SetLineWidth(width)
        a.SetLabelSize(lab_size if not line else 0)
        a.SetTitleSize(tit_size)
        a.SetTitleOffset(off)
        a.SetTitle(title)
        a.SetTitleColor(col)
        a.SetLabelColor(col)
        a.SetLabelFont(42)
        a.SetTitleFont(42)
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        self.Objects.append(a)
        return a

    def draw_y_axis(self, x, ymin, ymax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(x, x, ymin, ymax, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_x_axis(self, y, xmin, xmax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(xmin, xmax, y, y, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_line(self, x1, x2, y1, y2, color=1, width=1, style=1, name='li'):
        line = TCutG(name, 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw('same')
        self.Objects.append(line)
        return line

    def draw_tline(self, x1, x2, y1, y2, color=1, width=1, style=1):
        line = TLine(x1, y1, x2, y2)
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw()
        self.Objects.append(line)
        return line

    def draw_box(self, x1, y1, x2, y2, color=1, width=1, style=1, fillstyle=None, name='box', show=True):
        box = TCutG(name, 5, array([x1, x1, x2, x2, x1], 'd'), array([y1, y2, y2, y1, y1], 'd'))
        box.SetLineColor(color)
        box.SetFillColor(color)
        box.SetLineWidth(width)
        box.SetLineStyle(style)
        box.SetFillStyle(fillstyle) if fillstyle is not None else do_nothing()
        if show:
            box.Draw('same')
        self.Objects.append(box)
        return box

    def draw_vertical_line(self, x, ymin, ymax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(x, x, ymin, ymax, color, w, style, name) if not tline else self.draw_tline(x, x, ymin, ymax, color, w, style)

    def draw_horizontal_line(self, y, xmin, xmax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(xmin, xmax, y, y, color, w, style, name) if not tline else self.draw_tline(xmin, xmax, y, y, color, w, style)

    def draw_tlatex(self, x, y, text, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None):
        latex = TLatex(x, y, text)
        self.format_text(latex, name, align, color, size, angle, ndc, font)
        latex.Draw()
        self.Objects.append(latex)
        return latex

    def draw_arrow(self, x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        self.Objects.append(ar)

    def draw_tpad(self, name, tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None):
        margins = [.1, .1, .1, .1] if margins is None else margins
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
        p.SetFillColor(fill_col)
        p.SetMargin(*margins)
        do([p.SetLogx, p.SetLogy, p.SetLogz], [logx, logy, logz])
        do([p.SetGridx, p.SetGridy], [gridx, gridy])
        make_transparent(p) if transparent else do_nothing()
        p.Draw()
        p.cd()
        self.Objects.append(p)
        return p

    def draw_tpavetext(self, text, x1, x2, y1, y2, font=42, align=0, size=0, angle=0, margin=.05, color=1):
        p = TPaveText(x1, y1, x2, y2, 'ndc')
        p.SetFillColor(0)
        p.SetFillStyle(0)
        p.SetBorderSize(0)
        p.SetMargin(margin)
        t = p.AddText(text)
        self.format_text(t, 'pave', align, color, size, angle, ndc=True, font=font)
        p.Draw()
        self.Objects.append(p)
        return p

    def draw_preliminary(self, canvas=None, height=.06):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        return self.draw_tpavetext('#font[62]{RD42} Preliminary', c.GetLeftMargin(), .5, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=72, align=12, margin=0.04)

    def draw_irradiation(self, irr, canvas=None, height=.06, left=True):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        x1, x2 = (c.GetLeftMargin(), .5) if left else (.5, 1 - c.GetRightMargin())
        return self.draw_tpavetext('Irradiation: {}'.format(irr), x1, x2, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=42, align=12, margin=0.04)

    def draw_stats(self, fit, y2=None, width=.3, prec='5.1f', names=None):
        names = fit.Names if names is None else names
        c = get_last_canvas()
        tm = .98 - .05 - c.GetTopMargin() if y2 is None else y2
        rm = .98 - c.GetRightMargin()
        p = TPaveStats(rm - width, tm - .06 * (fit.NPars + 1), rm, tm, 'ndc')
        p.SetBorderSize(1)
        p.SetFillColor(0)
        p.SetFillStyle(0)
        latex = p.AddText('Fit Result')
        latex.SetTextFont(42)
        ls = p.GetListOfLines()
        ls.Add(self.draw_tlatex(0, 0, '#chi^{{2}} / ndf  = {chi2:{p}} / {ndf}'.format(ndf=fit.Ndf(), chi2=fit.Chi2(), p=prec), size=0, align=0, font=42))
        for i in xrange(fit.NPars):
            ls.Add(self.draw_tlatex(0, 0, '{n}  = {v:{p}} #pm {e:{p}}'.format(n=names[i], v=fit.Parameter(i), e=fit.ParError(i), p=prec), size=0, align=0, font=42))
        p.Draw()
        self.Objects.append(p)
        return p

    def draw_frame(self, pad, xmin, xmax, ymin, ymax, tit, div=None, y_cent=None):
        pad.cd()
        fr = pad.DrawFrame(xmin, ymin, xmax, ymax)
        pad.Modified()
        fr.GetYaxis().SetTitle(tit)
        do(fr.GetYaxis().CenterTitle, y_cent)
        fr.GetYaxis().SetNdivisions(div) if div is not None else do_nothing()
        self.format_frame(fr)
        self.Objects.append(fr)
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region SAVING

    def save_plots(self, savename, sub_dir=None, canvas=None, prnt=True, save=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        canvas = get_last_canvas() if canvas is None else canvas
        canvas.Modified()
        canvas.Update()
        if save:
            try:
                self.save_canvas(canvas, sub_dir=sub_dir, name=savename, print_names=prnt)
                self.Objects.append(canvas)
            except Exception as inst:
                warning('Error in save_canvas:\n{0}'.format(inst))

    def save_canvas(self, canvas, sub_dir=None, name=None, print_names=True):
        """ Saves the provided canvas into all the FileTypes. """
        sub_dir = self.SaveDir if sub_dir is None else sub_dir
        file_name = canvas.GetName() if name is None else name
        file_path = join(self.ResultsDir, sub_dir, '{typ}', file_name)
        set_root_output(False)
        for ext in self.FileTypes:
            ensure_dir(dirname(file_path.format(typ=ext)))
            canvas.SaveAs('{f}.{ext}'.format(f=file_path, ext=ext).format(typ=ext))
        if print_names:
            info('Saving plots: {nam}'.format(nam=file_name), prnt=self.Verbose)
        set_root_output(True)

    def save_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=.1, rm=.03, bm=.1, tm=None, draw_opt='', x_fac=None, y_fac=None, leg=None, logy=None, logx=None,
                   logz=None, canvas=None, grid=False, gridx=False, gridy=False, save=True, prnt=True, phi=None, theta=None):
        tm = (.1 if self.ActivateTitle else .03) if tm is None else tm
        x = self.Res if x_fac is None else int(x_fac * self.Res)
        y = self.Res if y_fac is None else int(y_fac * self.Res)
        h = histo
        set_root_output(show)
        c = TCanvas('c_{0}'.format(h.GetName()), h.GetTitle().split(';')[0], x, y) if canvas is None else canvas
        c.SetMargin(lm, rm, bm, tm)
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx or grid, gridy or grid])
        do([c.SetPhi, c.SetTheta], [phi, theta])
        h.Draw(draw_opt)
        if leg is not None:
            for i in [leg] if type(leg) is not list else leg:
                i.Draw()
        self.save_plots(save_name, sub_dir=sub_dir, prnt=prnt, save=save)
        set_root_output(True)
        self.Objects.append([c, h, leg] if leg is not None else [c, h])
        return c

    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region FORMATTING
    @staticmethod
    def format_text(t, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None):
        t.SetName(name)
        t.SetTextAlign(align)
        t.SetTextColor(color)
        t.SetTextSize(size)
        do(t.SetTextAngle, angle)
        do(t.SetTextFont, font)
        do(t.SetNDC, ndc)
        return t

    @staticmethod
    def format_histo(histo, name=None, title=None, x_tit='', y_tit='', z_tit='', marker=20, color=None, markersize=None, x_off=None, y_off=None, z_off=None, lw=1,
                     fill_color=None, fill_style=None, stats=True, tit_size=None, lab_size=None, l_off_y=None, l_off_x=None, draw_first=False, x_range=None, y_range=None, z_range=None,
                     do_marker=True, style=None, ndivx=None, ndivy=None, ncont=None, tick_size=None, t_ax_off=None, center_y=False, center_x=False, yax_col=None):
        h = histo
        if draw_first:
            set_root_output(False)
            h.Draw('a')
            set_root_output(True)
        do(h.SetTitle, title)
        do(h.SetName, name)
        try:
            h.SetStats(stats)
        except AttributeError or ReferenceError:
            pass
        # markers
        try:
            if do_marker:
                do(h.SetMarkerStyle, marker)
                do(h.SetMarkerColor, color)
                do(h.SetMarkerSize, markersize)
        except AttributeError or ReferenceError:
            pass
        # lines/fill
        try:
            h.SetLineColor(color) if color is not None else h.SetLineColor(h.GetLineColor())
            h.SetLineWidth(lw)
            h.SetFillColor(fill_color) if fill_color is not None else do_nothing()
            h.SetFillStyle(fill_style) if fill_style is not None else do_nothing()
            h.SetFillStyle(style) if style is not None else do_nothing()
            h.SetContour(ncont) if ncont is not None else do_nothing()
        except AttributeError or ReferenceError:
            pass
        # axis titles
        try:
            x_axis = h.GetXaxis()
            if x_axis:
                x_axis.SetTitle(x_tit) if x_tit else h.GetXaxis().GetTitle()
                x_axis.SetTitleOffset(x_off) if x_off is not None else do_nothing()
                do(x_axis.SetTitleSize, tit_size)
                do(x_axis.SetLabelSize, lab_size)
                x_axis.SetRangeUser(x_range[0], x_range[1]) if x_range is not None else do_nothing()
                x_axis.SetNdivisions(ndivx) if ndivx is not None else do_nothing()
                do(x_axis.SetLabelOffset, l_off_x)
                do(x_axis.SetTickSize, tick_size)
                x_axis.CenterTitle(center_x)
            y_axis = h.GetYaxis()
            if y_axis:
                y_axis.SetTitle(y_tit) if y_tit else y_axis.GetTitle()
                y_axis.SetTitleOffset(y_off) if y_off is not None else do_nothing()
                do(y_axis.SetTitleSize, tit_size)
                do(y_axis.SetLabelSize, lab_size)
                do(y_axis.SetLabelOffset, l_off_y)
                y_axis.SetRangeUser(y_range[0], y_range[1]) if y_range is not None else do_nothing()
                do(y_axis.SetNdivisions, ndivy)
                y_axis.CenterTitle(center_y)
                do(y_axis.SetTitleColor, yax_col)
                do(y_axis.SetLabelColor, yax_col)
                do(y_axis.SetAxisColor, yax_col)
            z_axis = h.GetZaxis()
            if z_axis:
                z_axis.SetTitle(z_tit) if z_tit else h.GetZaxis().GetTitle()
                z_axis.SetTitleOffset(z_off) if z_off is not None else do_nothing()
                do(z_axis.SetTitleSize, tit_size)
                do(z_axis.SetLabelSize, lab_size)
                z_axis.SetRangeUser(z_range[0], z_range[1]) if z_range is not None else do_nothing()
        except AttributeError or ReferenceError:
            pass
        set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()

    @staticmethod
    def format_pie(pie, h=None, r=None, text_size=None, angle3d=None, angle_off=None, label_format=None):
        do([pie.SetHeight, pie.SetRadius], [h, r])
        do(pie.SetTextSize, text_size)
        do(pie.SetAngle3D, angle3d)
        do(pie.SetLabelFormat, label_format)
        do(pie.SetAngularOffset, angle_off)

    def format_statbox(self, x=.95, y=None, w=.16, entries=3, only_fit=False, fit=False, only_entries=False, opt=None, form=None):
        y = .88 if self.ActivateTitle else .95 if y is None else y
        if only_fit or fit:
            gStyle.SetOptStat(1111 if fit else 0011)
            gStyle.SetOptFit(1)
        if only_entries:
            gStyle.SetOptStat(1000000010 if not only_fit else 1000000011)
            entries = 6 if entries == 3 else entries
        gStyle.SetOptStat(opt) if opt is not None else do_nothing()
        gStyle.SetFitFormat(form) if form is not None else do_nothing()
        gStyle.SetStatX(x)
        gStyle.SetStatY(y)
        gStyle.SetStatW(w)
        gStyle.SetStatH(.02 * entries)

    @staticmethod
    def format_frame(frame):
        fr = frame
        fr.GetYaxis().SetTitleSize(.06)
        fr.GetYaxis().SetTitleOffset(.6)
        fr.GetYaxis().SetLabelSize(.06)
        fr.SetTitleSize(.05)
        fr.GetXaxis().SetTickLength(0)
        fr.GetXaxis().SetLabelOffset(99)
        fr.SetLineColor(0)
        fr.GetXaxis().SetTimeDisplay(1)
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region MAKING
    @staticmethod
    def make_tgrapherrors(name, title, color=1, marker=20, marker_size=1, width=1, asym_err=False, style=1, x=None, y=None, ex=None, ey=None):
        x = list(x) if type(x) == ndarray else x
        if x is None:
            gr = TGraphErrors() if not asym_err else TGraphAsymmErrors()
        else:
            gr = TGraphErrors(*make_graph_args(x, y, ex, ey)) if not asym_err else TGraphAsymmErrors(*make_graph_args(x, y, ex, ey))
        gr.SetTitle(title)
        gr.SetName(name)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        gr.SetMarkerSize(marker_size)
        gr.SetLineWidth(width)
        gr.SetLineStyle(style)
        return gr

    def make_legend(self, x1=.65, y2=.88, nentries=2, scale=1, name='l', y1=None, clean=False, margin=.25, x2=None):
        x2 = .95 if x2 is None else x2
        y1 = y2 - nentries * .05 * scale if y1 is None else y1
        legend = TLegend(x1, y1, x2, y2)
        legend.SetName(name)
        legend.SetTextFont(42)
        legend.SetTextSize(0.03 * scale)
        legend.SetMargin(margin)
        if clean:
            legend.SetLineWidth(2)
            legend.SetBorderSize(0)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetTextAlign(12)
        self.Objects.append(legend)
        return legend

    def make_canvas(self, name='c', title='c', x=1., y=1., show=True, logx=None, logy=None, logz=None, gridx=None, gridy=None, transp=None):
        set_root_output(show)
        c = TCanvas(name, title, int(x * self.Res), int(y * self.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        self.Objects.append(c)
        return c

    def make_graph_from_profile(self, p):
        x_range = [i for i in xrange(p.GetNbinsX()) if p.GetBinContent(i)]
        x = [make_ufloat([p.GetBinCenter(i), p.GetBinWidth(i) / 2]) for i in x_range]
        y = [make_ufloat([p.GetBinContent(i), p.GetBinError(i)]) for i in x_range]
        return self.make_tgrapherrors('g{n}'.format(n=p.GetName()[1:]), p.GetTitle(), x=x, y=y)
    # endregion
    # ----------------------------------------


def create_colorlist():
    col_names = [kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, kCyan, kTeal]
    colors = []
    for color in col_names:
        colors.append(color + (1 if color != 632 else -7))
    for color in col_names:
        colors.append(color + (3 if color != 800 else 9))
    return colors


def make_graph_args(x, y, ex=None, ey=None):
    if len(list(x)) != len(list(y)):
        warning('Arrays have different size!')
        return []
    lx = len(x)
    if type(x[0]) is Variable:
        return [lx, array([v.n for v in x], 'd'), array([v.n for v in y], 'd'), array([v.s for v in x], 'd'), array([v.s for v in y], 'd')]
    return [lx, array(x, 'd'), array(y, 'd'), array(ex, 'd') if ex is not None else zeros(lx), array(ey, 'd') if ey is not None else zeros(lx)]


def make_transparent(pad):
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        warning('There is no canvas is in the list...')


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


def set_time_axis(histo, form='%H:%M', off=0):
    histo.GetXaxis().SetTimeFormat(form)
    histo.GetXaxis().SetTimeOffset(-off - 3600 if off else 0)
    histo.GetXaxis().SetTimeDisplay(1)


def set_palette(pal):
    gStyle.SetPalette(pal)


def make_ufloat(tup):
    if type(tup) is Variable:
        return tup
    return ufloat(tup[0], tup[1])


def set_statbox(x=.95, y=.88, w=.16, entries=3, only_fit=False, fit=False, only_entries=False, opt=None, form=None):
    if only_fit or fit:
        gStyle.SetOptStat(1111 if fit else 0011)
        gStyle.SetOptFit(1)
    if only_entries:
        gStyle.SetOptStat(1000000010 if not only_fit else 1000000011)
        entries = 6 if entries == 3 else entries
    gStyle.SetOptStat(opt) if opt is not None else do_nothing()
    gStyle.SetFitFormat(form) if form is not None else do_nothing()
    gStyle.SetStatX(x)
    gStyle.SetStatY(y)
    gStyle.SetStatW(w)
    gStyle.SetStatH(.02 * entries)


def set_z_range(zmin, zmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetZaxis().SetRangeUser(zmin, zmax)


def load_resolution():
    global resolution
    if resolution is None:
        try:
            resolution = round_down_to(get_monitors()[0].height, 500)
        except Exception as err:
            warning(err)
            return 1000
    return resolution
