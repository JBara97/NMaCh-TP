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

from Data_processing import df_TEs, df_interactions

df_TEs_sorted = df_TEs.sort_values(by=["Threshold (°C) – Est."])
tps = df_TEs_sorted[['Threshold (°C) – Min.','Threshold (°C) – Est.','Threshold (°C) – Max.']].to_numpy()

# Pre-making the default Basemap
m = Basemap(projection='kav7',lon_0=-30,resolution='c',)

# Pre-making the default graph from df_interactions
G = nx.from_pandas_edgelist(df_interactions,source='Source',target='Target',edge_attr=['min', 'max'],create_using=nx.DiGraph)

vmin = -1
vmax = 9 
node_cmap = 'rocket'
node_size=[1600/np.log(df_TEs.loc[n]['Collapse (years)']) for n in G.nodes]

earc_dic = {e:0 for e in G.edges}
earc_dic[('GIS','AMOC')]=0.2
earc_dic[('AMOC','GIS')]=0.2
earc_dic[('AMOC','WAIS')]=-0.5
earc_dic[('AMOC','ENSO')]=0.3
earc_dic[('WAIS','AMOC')]=0.3
earc_dic[('ASI','AMOC')]=0.15
earc_dic[('AMOC','ASI')]=0.15
earc_dic[('ASI','GIS')]=0.1
earc_dic[('PERM','AMOC')]=-0.15
earc_dic[('ASI','PERM')]=-0.1
earc_dic[('PERM','ASI')]=-0.1
earc_dic[('WAIS','GIS')]=-0.9
earc_dic[('GIS','WAIS')]=1.4

dest_edge_color = 'indianred'
stab_edge_color = 'royalblue'
uncl_edge_color = 'gray'

pos_lonlat = {'GIS':[-43,75], 'AMOC':[-30,40],'ENSO':[-115,-8.8],
                  'WAIS':[-90,-75],'ASI':[0,68],'PERM':[70,60],
                  'WAM':[0,15],'REEF':[-80,12],'AMAZ':[-62.2,-3.5]}
pos_xy = {k:m(v[0],v[1]) for k,v in pos_lonlat.items()}

xmax,ymax = (34667130.92439352,20015077.371242613)

T_lim_dic = df_TEs['Threshold (°C) – Est.'].to_dict()

node_size=[1600/np.log(df_TEs.loc[n]['Collapse (years)']) for n in G.nodes]
node_labels = {'GIS':'Greenland ice sheet (GIS)', 'AMOC':'AMOC','ENSO':'ENSO',
              'WAIS':'West Antarctic ice sheet (WAIS)','ASI':'Arctic sea ice (ASI)','PERM':'Permafrost collapse (PERM)',
              'WAM':'West African monsoon (WAM)','REEF':'Coral reefs (REEF)','AMAZ':'Amazon dieback (AMAZ)'}
node_label_offsets = {'GIS': np.array([0,0.035]), 
                      'AMOC':np.array([-0.025,0.01]),
                      'ENSO':np.array([-.025,0]),
                      'WAIS':np.array([0,-0.035]),
                      'ASI':np.array([0.025,0.025]),#np.array([0,0.05]),
                      'PERM':np.array([0.025,0]),
                      'WAM':np.array([0,0.05]),
                      'REEF':np.array([0,0.05]),
                      'AMAZ':np.array([0,-0.042])}
node_label_pos = {k:np.array(v) +
                  node_label_offsets[k]*np.array([xmax,ymax])
                  for k,v in pos_xy.items()}
node_label_ha = {k:'center' for k in G.nodes}
node_label_ha['AMOC'] = 'right'
node_label_ha['ENSO'] = 'right'
node_label_ha['WAM'] = 'left'
node_label_ha['GIS'] = 'center'
node_label_ha['ASI'] = 'left'
node_label_ha['PERM'] = 'left'

node_label_va = {k:'bottom' if v[1]>0 else 'top' if v[1]<0 else 'center' 
                 for k,v in node_label_offsets.items()}

bbox_node=dict(facecolor='white', 
               edgecolor='k', 
               boxstyle='round,pad=0.2')
bbox_white=dict(facecolor='white', 
               edgecolor='white', 
               boxstyle='round,pad=0.2')

def network_panel(fig,ax,
                  node_size=node_size,
                  node_cmap='rocket',
                  vmin=-1,
                  vmax=9,
                  m = m,
                  G = G,
                  earc_dic = earc_dic,
                  pos_xy = pos_xy, 
                  T_lim_dic = T_lim_dic,
                  node_labels = node_labels,
                  node_label_offsets=node_label_offsets,
                  node_label_ha = node_label_ha,
                  node_label_va = node_label_va,
                  bbox_node = bbox_node,
                  dest_edge_color = 'indianred',
                  stab_edge_color = 'royalblue',
                  uncl_edge_color = 'gray'
                 ):  
    def make_legend_arrow(legend, orig_handle,
                          xdescent, ydescent,
                          width, height, fontsize):
        p = mpl.patches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
        return p
    
    m.fillcontinents(color='tan',lake_color='tan',ax=ax,zorder=1);
    m.drawcoastlines(color='k',zorder=0,ax=ax)    
    
    nx.draw_networkx(G,ax=ax,pos=pos_xy,
                     node_size=node_size,with_labels=False,edgelist=[],
                     node_color=[T_lim_dic[n] for n in G.nodes()],
                     vmin=vmin,vmax=vmax,cmap=node_cmap)
    for e in G.edges:
        vals = G.get_edge_data(e[0],e[1])
        if vals['min'] > 0: # Destabilising
            c = 'indianred'
            w = vals['max']
        elif vals['max'] < 0: # Stabilising
            c = 'royalblue'
            w= np.abs(vals['min'])
        else: # Unclear
            c = 'gray' 
            w = vals['max']
        nx.draw_networkx_edges(G,pos_xy,ax=ax,edgelist=[e], 
                               edge_color=c, width=np.log(w)*1.8,
                               node_size=[n+200 for n in node_size],
                               connectionstyle="arc3,rad=rrr".replace('rrr',str(earc_dic[e])),
                               arrowsize=20
                              );
    for n,(x,y) in pos_xy.items():
        ax.annotate(node_labels[n], (x,y), node_label_pos[n],
                    ha=node_label_ha[n],
                    va=node_label_va[n],
                    bbox=bbox_node)
    
    custom_lines = [FancyArrow(0,0,5,0, color=dest_edge_color,
                               label='Destabilising',shape='right',
                               length_includes_head=True),
                    FancyArrow(0,0,5,0, color=stab_edge_color,
                               label='Stabilising',shape='right',
                               length_includes_head=True),
                    FancyArrow(0,0,5,0, color=uncl_edge_color,
                               label='Unclear',shape='right',
                               length_includes_head=True)]
    ax.legend(handles=custom_lines,
              handler_map={mpl.patches.FancyArrow : 
                           mpl.legend_handler.HandlerPatch(patch_func=make_legend_arrow)},
              loc='lower left',title='Interaction type',edgecolor='k',framealpha=1)
    return ax


def threshold_inset(fig,ax,
                    tps=tps, 
                    TEs_labels = df_TEs_sorted.index.values,
                    bbox = bbox_white,
                    xlim=[0,8.3],
                    ylim=[-.6,8.6],
                   ):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.barh(np.arange(9),[z-x for x,y,z in tps],
            left=[x for x,y,z in tps]) # horizontal
    ax.set_yticks(np.arange(9),TEs_labels)
    for i,v in enumerate(tps):
        x,y,z=v
        im=ax.imshow(np.linspace(x,z).reshape((1,50)), extent=[x,z,i-0.4,i+0.4], aspect="auto", 
                     zorder=1,vmin=vmin,vmax=vmax,cmap=node_cmap,origin='lower')
        ax.plot([y]*2,[i-0.4,i+0.4],c='k',ls=':')
    
    ax.set_xticks([])
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',pad=0)
    cb.ax.set_xlim(xlim)
    cb.ax.set_xticks(np.arange(9))
    cb.ax.tick_params(axis='x',direction='in')
    cb.ax.set_xlabel('Temperature (°C)',bbox=bbox)
    
    ax.yaxis.tick_right()
    ax.tick_params(axis="y",direction="in")
    ax.patch.set_alpha(0.7)
    plt.setp(ax.get_yticklabels(), 
             bbox=bbox)
    plt.setp(cb.ax.get_xticklabels(), 
             bbox=bbox)


def threshold_panel(fig,ax,
                    cmap='rocket',
                    vmin=-1,
                    vmax=9,
                    ylim=[0,9],
                    tps=tps, 
                    TEs_labels = df_TEs_sorted.index.values,):
    ax.bar(np.arange(9),[z-x for x,y,z in tps],bottom=[x for x,y,z in tps])
    ax.set_ylim(ylim)
    
    for i,v in enumerate(tps):
        x,y,z=v
        ax.imshow(np.linspace(x,z).reshape((50,1)), extent=[i-0.4,i+0.4,x,z], aspect="auto", 
                  zorder=1,vmin=vmin,vmax=vmax,cmap=cmap,origin='lower')
        ax.plot([i-0.4,i+0.4],[y]*2,c='k')
    
    ax.set_xlim(-1,9)
    ax.set_xticks(np.arange(9),TEs_labels)
    ax.tick_params("x", rotation=45)


def stylisedOS_panel(fig,ax,
                     stylized_par=(1.2,3.5,2,80,20,100),
                     temp_names = [r'$\Delta T_0$', r'$\Delta T_p$', r'$\Delta T_s$'],
                     dur_names = [r'$t_u$', r'$t_p$', r'$t_d$'],
                     phase_names = ['Ramp up', 'Peak', 'Ramp down', 'Stabilisation']
                    ):
    ylim=ax.get_ylim()
    xlim=ax.get_xlim()
    xdata = np.arange(300)
    dashed_ends = [0,stylized_par[3],np.sum(stylized_par[3:])]#[0,popt_of_m[3],np.sum(popt_of_m[3:])]
    for i in range(3):
        ax.plot([-20,dashed_ends[i]],[stylized_par[i]]*2,c='k',ls='--',lw=0.5)
        ax.annotate('',
                    (np.sum(stylized_par[3:3+i]),4), 
                    (np.sum(stylized_par[3:4+i]),4),
                    arrowprops=dict(facecolor='black',arrowstyle='<|-|>,head_width=0.4,head_length=0.5')
                   )
        ax.text(np.sum(stylized_par[3:3+i])+stylized_par[3+i]/2, 3.7,
                dur_names[i],va='bottom',ha='center')
        ax.text(np.sum(stylized_par[3:3+i])+stylized_par[3+i]/2, 4.5,
                phase_names[i],va='bottom',ha='center')
    ax.text(np.mean([xlim[1],np.sum(stylized_par[3:])]), 4.5,
                phase_names[-1],va='bottom',ha='center')
    ax.annotate('',
                (np.sum(stylized_par[3:]),4), 
                (xlim[1],4),
                arrowprops=dict(facecolor='black',arrowstyle='<-')
               )
    
    ax.plot(xdata,OS_model(xdata,*stylized_par),color='#607dab',lw=2)#*np.median(popts,axis=0)))
    
    ax.fill_betweenx(ylim,stylized_par[3],np.sum(stylized_par[3:5]),color='silver',zorder=0)
    ax.fill_betweenx(ylim,np.sum(stylized_par[3:]),xlim[1],color='silver',zorder=0)
    
    ax.set_xlabel("Model time (yr)")
    ax.set_ylabel('Temperature (°C)')
    
    ax.set_yticks(stylized_par[:3],temp_names)