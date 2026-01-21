### Utility functions for plotting
import matplotlib.pyplot as plt
import colorsys
import matplotlib.colors as mc
import numpy as np

from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

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

def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` Axes.

    This function creates a RadarAxes projection and registers it.

    Function from Matplotlib documentation: https://matplotlib.org/stable/gallery/specialty_plots/radar_chart.html

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding Axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta