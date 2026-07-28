import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow,FancyArrowPatch
from matplotlib.gridspec import GridSpec
import networkx as nx
import pandas as pd
from mpl_toolkits.basemap import Basemap
import numpy as np
import seaborn as sns
from util.overshoots import OS_model
from util.plot import lighten_color
from matplotlib.colors import LinearSegmentedColormap
import colorcet as cc

colors_g = mpl.colormaps['Greens'](np.linspace(0, 1, 12))
colors_h = mpl.colormaps['cet_fire_r'](np.linspace(0, 1, 12))
c_offset = 1
c_list = [colors_g[7]] + [c_h for c_h in colors_h[c_offset:c_offset+9]]
cmap = LinearSegmentedColormap.from_list('my_cmap', c_list, N=100).with_extremes(bad='white')

def risk_heatmap(x,y,risk,
                 cmap=cmap,label='Tipping risk (%)',
                 clevels=[1,5]+[i for i in range(10,100,10)],ccols=['k'],
                 vmin=0,vmax=100,ax=None):
    if ax is None:
        ax=plt.gca()
    X,Y=np.meshgrid(x,y,indexing='xy')
    CS = ax.contour(X,Y, risk,clevels,
                    corner_mask=True,                    
                    colors=ccols,linewidths=0.5,zorder=1,ls='-')
    ax.clabel(CS)
    return ax.imshow(risk,cmap=cmap,#cmap='Reds',
                 extent=(x[0] - np.diff(x)[0]/2, 
                         x[-1]+ np.diff(x)[0]/2,
                         y[0] - np.diff(y)[0]/2, 
                         y[-1]+ np.diff(y)[0]/2),
                 origin='lower',
                 vmin=vmin,vmax=vmax,zorder=0,
                 aspect='auto')