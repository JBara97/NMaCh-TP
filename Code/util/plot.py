### Utility functions for plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import colorsys
import matplotlib.colors as mc
import numpy as np

from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
import colorcet

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple. 
    If amount >> 1 then an error might be raised as the corresponding hex string would not be valid

    Parameters
    ---
    color : color string, hex string, or RGB tuple. 
        The color to be modified
    amount : float
        The amount to lighten the colour if amount is in [0,1] otherwise if amount > 1 then the colour is darkened instead

    Returns
    ---
    out : RGB colour
    
    Examples
    ---
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def TE_3x3_fig(xlabel='Year',ylabel='Degraded prob. (%)',figsize=(12,9),sharex=True,sharey=False, 
               TE_ordered=['AMOC', 'WAM', 'REEF',
                           'ENSO', 'AMAZ','PERM',
                           'GIS', 'WAIS', 'ASI'], 
               subfig_lab_loc=(0.95,0.05), ha='right',va='bottom',
               fontsize=12, 
               bbox=dict(fc='white', edgecolor='0.7',boxstyle='round')):
    """
    Produces a 3 x 3 figure where each subfigure is labelled with the tipping element name. 
    Moreover the xlabel(s) and ylabel(s) will be displayed

    Parameters (optional)
    ---
    xlabel, ylabel : str
    figsize : 2-tupple
    sharex,sharey : Bool
    TE_ordered : list or array
        Ordered list of tipping elements, this will form the label for each subfigure
    subfig_lab_loc : tupple or list of tupples
        If tupple must be a 2-tupple which gives the (x,y)-position for all subfigure labels
        Else, must be a list of such 2-tupple positions
    ha,va : str or list of str
        If str gives the shared text alignment for horizontal (ha) or vertical (va).
        If list of str must be of 9 such strings.
    fontsize : number
    bbox : dictionary for bbox attributes for the subfigure labels

    Returns
    ---
    fig : matplotlib.figure.Figure
    axs : numpy.ndarray
        3 x 3 array where each element is a matplotlib.axes._axes.Axes
    """
    fig,axs = plt.subplots(figsize=figsize,sharex=sharex,sharey=sharey,ncols=3,nrows=3)
    fig.subplots_adjust(hspace=0.05)
    letters = ['(%s)'%l for l in 'abcdefghi']
    if len(subfig_lab_loc) == 2:
        subfig_lab_loc = [subfig_lab_loc]*9
    if type(ha)==str:
        ha=[ha]*9
    if type(va)==str:
        va=[va]*9
    for i in range(9):
        axs[i//3,i%3].grid(visible=True,lw=0.5,axis='y')
        axs[i//3,i%3].annotate(letters[i] + ' '+ TE_ordered[i],
                               subfig_lab_loc[i],
                               xycoords='axes fraction',
                               ha=ha[i], va=va[i],
                               fontsize=fontsize,
                               bbox=bbox)
    for i in range(3):
        axs[-1,i].set_xlabel(xlabel)
        axs[i,0].set_ylabel(ylabel)
    return fig,axs

def ribbon_plot(ax,x,y,quantiles=[0.05],axis=0,color='k',shade=None,alpha=None,
                edgecolor='none',label=None,zorder=None):
    """
    Plots ribbon plot(s) on axis ax for x-values x and data y.
    For each quantile in quantiles:
        a ribbon plot will be produced whose lower value is the quantile and upper value is 1-quantile, 
        i.e. symmetric quantiles around 0.5

    Parameters (required)
    ---
    ax : matplotlib.axes._axes.Axes
        Axis to draw in
    x : 1D array-like
        x-values. Must have the same length as one of the axes of y (specifically not the one in "axis")
    y : numpy.ndarray
        y-values. Needs to have at least 2 dimensions:
            one to plot over which has the same length as "x"
            the other(s) are the dimensions over which the quantiles are taken
    
    Parameters (optional)
    ---
    axis : int or tupple of int
        the axis or axes over which the quantiles will be taken
    quantiles : tupple/list of floats in [0,1]
        List of lower quantiles to plot. The upper values are calculated automatically by 1-quantile
    color : str
    shade, alpha : float
        list of shades/alpha-values. Each value corresponds to a ribbon plot specified by quantiles. 
        Must have the same length as quantiles
    edgecolor : str
        shared edgecolor for all ribbon plots
    """
    if alpha is None:
        alpha = [1]*len(quantiles)
    if shade is None:
        shade = quantiles
    if zorder is None:
        zorder = np.arange(len(quantiles))
    if label is None:
        label = ['($P_{%s}, P_{%s}$)'%(int(q*100),int(100-q*100)) for q in quantiles]
    for i in range(len(quantiles)):
        ax.fill_between(x, np.quantile(y,quantiles[i], axis=axis), 
                        np.quantile(y,1-quantiles[i],axis=axis)  , 
                        color=lighten_color(color,shade[i]),
                        alpha=alpha[i], edgecolor=edgecolor,label=label[i],zorder=zorder[i])
    return

class Radar(object):
    def __init__(self, fig, titles, rect=None, wedge_width=None,rotate_labels=True, with_rlabels=True, loc_offset=-2,
                 rticks=None, ylims = None):
        if rect is None:
            rect = [0.05, 0.05, 0.95, 0.95]

        self.n = len(titles)
        self.angles = [a if a <=360. else a - 360. for a in np.arange(90, 90+360, 360.0/self.n)] # in degrees
        self.axes = [fig.add_axes(rect, projection="polar", label="axes%d" % i) 
                         for i in range(self.n)]
        self.rticks = rticks
        if self.rticks is None:
            self.rticks = [ax.get_yticks() for ax in self.axes]
        self.ylims = ylims
        if self.ylims is None:
            self.ylims = [ax.get_ylim() for ax in self.axes]

        self.ax = self.axes[0]
        self.ax.set_thetagrids(self.angles, labels=titles,)# fontsize=12, weight="bold")
        self.ax.tick_params(axis='x', which='major', pad=16)
        self.ax.yaxis.grid(False)
        self.wedge_width = wedge_width
        if self.wedge_width is None:
            self.wedge_width = 360/self.n # width of each wedge in degrees

        if rotate_labels:
            labels = []
            for label, angle in zip(self.ax.get_xticklabels(),self.angles):
                x,y = label.get_position()
                lab = self.ax.text(x,y, label.get_text(), transform=label.get_transform(),
                              ha='center', va='center')
                lab.set_rotation(angle-90)
                labels.append(lab)
            self.ax.set_xticklabels([])

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)    
            ax.yaxis.grid(False)

        self.draw_rticks(rticks,loc_offset=loc_offset,width=self.wedge_width,with_label=with_rlabels)
        self.draw_wedgeticks(ylims=ylims,width=self.wedge_width)
    
    ### ------ Radar methods --------

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)

    def bar(self,values,angles=None,*args,**kw):
        if angles is None:
            angles = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        for ax, angle, value in zip(self.axes,angles, values):
            ax.bar(angle,value, *args, **kw)

    def draw_rticks(self,rticks, loc_offset=0, width=None, with_label=True):
        if width == None:
            width=self.wedge_width
        else:
            self.wedge_width = width
        for ax, angle, rtick in zip(self.axes, self.angles, rticks):
            ax.set_rgrids(rtick, labels = '', angle=angle, fontsize=12)
            ax.spines["polar"].set_visible(False)
            ax.set_rticks(rtick)

            ax.set_thetagrids([])
            theta1 = np.deg2rad(angle-width/2)
            theta2 = theta1 + np.deg2rad(width)
            for r in ax.get_yticks():
                ax.plot(np.linspace(theta1,theta2), [r]*50, color='k',lw = 0.7,zorder=0)
            if with_label:
                for rt in rtick:
                    ax.text(theta1+np.deg2rad(loc_offset),rt,rt,ha='center',va='center')
        self.rticks = rticks

    def draw_wedgeticks(self,ylims=None,width=None):
        if ylims is None:
            ylims = [ax.get_ylim() for ax in self.axes]
        if width is None:
            width=self.wedge_width
        for ax, angle, ylim in zip(self.axes, self.angles, ylims):
            # ylim=ax.get_ylim()
            ax.plot([np.deg2rad(angle - width/2)]*2, ylim, color='k',lw = 1,zorder=0)
            ax.plot([np.deg2rad(angle + width/2)]*2, ylim, color='k',lw = 1,zorder=0)
            ax.set_ylim(ylim)

    def wedge_stacked_bar(self, ax, angle, values, 
                          width=None, c_list=None, 
                          degrees=True, bottom = None, 
                          edgecolor='none', alpha=0.7,
                          **kwargs):
        if width is None:
            width=np.deg2rad(self.wedge_width)
        if degrees:
            angle=np.deg2rad(angle)
        if c_list is None:
            c_list = ['g'] + [c for c in mpl.colormaps['cet_fire_r'](np.linspace(0, 1, len(values)+4))[1:len(values)]]
        
        if bottom is None:
            bottom = ax.get_ylim()[0]
        for v,c in zip(values,c_list[:len(values)]):
            ax.bar(angle, v-bottom, 
                   color=c, 
                   edgecolor=edgecolor,
                   alpha=alpha,
                   width = width,
                   bottom=bottom, **kwargs)
            bottom=max(bottom,v)