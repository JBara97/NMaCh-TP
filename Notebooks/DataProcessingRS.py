import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib

# Defining functions to convert strings into relevant numbers
def tword_to_number(tword):
    # NB: weeks = 2 weeks, months ~6 months (half a year), years ~ 5 years, decades ~ 5 decades, 
        # centuries ~ 5 centuries, unknown ~ 650 years
    # Here treat 1y = 365 days (not bothering too much about leap years)
    word_number = dict(Weeks=2*7, Months=365/2, Years=5*365, Decades=50*365, 
                       Centuries=500*365, Unknown=600*365)
    return word_number[tword]

# Turning a tstring like 'Weeks, Months, Years' into a list of day-equivalents eg [2*7, 365.25/2, 5*365.25]
def tstring_to_numberslist(tstr):
    twords = []
    commaix = tstr.find(',')
    while commaix>0:
        tword=tstr[: commaix]
        tstr=tstr[commaix+2 :]
        commaix = tstr.find(',')
        twords.append(tword)
    twords.append(tstr)
    return list(map(tword_to_number,twords))

# Turning a tstring like 'Weeks, Months, Years' into the exp(mean(log( ...) )) 
def tstring_to_number(tstr):
    return np.exp(np.mean(np.log(tstring_to_numberslist(tstr))))

# Turning a reversability string eg 'Irreversible, Unknown' into the corresponding days-equivalent
def revstring_to_number(revstr, collapsestr):
    twords = []
    numbers = []
    commaix = revstr.find(',')
    while commaix>0:
        tword=revstr[: commaix]
        revstr=revstr[commaix+2 :]
        commaix = revstr.find(',')
        twords.append(tword)
    twords.append(revstr)
    
    for t in twords:
        if 'Readily reversible' in t:
            # Readily reversible => smallest collapse time scale available
            numbers.append(tstring_to_numberslist(collapsestr)[0])
        elif 'Irreversible' in t:
            # Irreversible => 7 centuries
            numbers.append(700*365.25)
        elif 'Hysteretic' in t:
            # Hysteretic => collapse time scale
            numbers.append(tstring_to_number(collapsestr))
        elif 'Unknown' in t:
            # Unknown => 6 centuries
            numbers.append(600*365.25)
        else: # If not in the above then use a -1 which might raise the error from a log(-1)
            numbers.append(-1)
    return np.exp(np.mean(np.log(numbers)),dtype=float)

# Splitting up a str of a list into a list of strs
def strlist_to_liststr(strlist):
    if type(strlist) == str:
        return strlist.split(', ')
    else:
        return []

### Reading & processing the Regime Shifts Database (RSD) as used by Rocha et al. 2018, i.e. the one accessed on the 25th Oct 2016
## Reading the database
fpath_old = '../Data/161025_RegimeShiftsDatabase.csv'
RSD = pd.read_csv(fpath_old,encoding="ISO-8859-1") # there are 32 regime shifts in here
## Reading the causal loop diagram (CLD) data (which is an edgelist) produced by Rocha et al 2018 
fpath = '../Data/RS_CLD_2018.csv'
df_edgelist = pd.read_csv(fpath, sep=';') # it's an edgelist!! containing 30 regime shift nodes
df_edgelist = df_edgelist[['Tail','Head','Polarity','Biome','Regime Shift']]
df_edgelist['Regime Shift'] = df_edgelist['Regime Shift'].str.lower() # lower casing regime shift names for semantic cohesion later
df_edgelist = df_edgelist[df_edgelist['Polarity'].isin(['-1','1'])] # following the methods of Rocha et al 2018

## Add desertification as is done in the Rocha et al 2018 SI code
desert_dict = {'Regime Shift Name':['Desertification'], 
               'Drivers - Key direct drivers':["Vegetation conversion and habitat fragmentation, Harvest and resource consumption, Infrastructure development, Soil erosion & land degradation, Environmental shocks (eg floods), Global climate change"],
               'Drivers - Land use': ["Small-scale subsistence crop cultivation, Large-scale commercial crop cultivation, Extensive livestock production (natural rangelands)"],
               'Impacts - Ecosystem type': ['Drylands & deserts (below ~500mm rainfall/year), Grasslands'],
               'Impacts - Key Ecosystem Processes': ["Soil formation, Primary production, Water cycling"],
               'Impacts - Biodiversity': ['Biodiversity'],
               'Impacts - Provisioning services':["Freshwater, Food Crops, Livestock, Wild animal and plant products, Woodfuel"],
               'Impacts - Regulating services':["Air quality regulation, Climate regulation, Regulation of soil erosion"],
               'Impacts - Cultural services':["Knowledge and educational values"],
               'Impacts - Human Well-being': ["Food and nutrition, Livelihoods and economic activity, Cultural, aesthetic and recreational values, Social conflict"],
               'Key Attributes - Typical spatial scale': ["Local/landscape (e.g. lake, catchment, community)"],
               'Key Attributes - Typical time scale': ["Years, Decades"],
               'Key Attributes - Reversibility': ["Hysteretic"],
               'Key Attributes - Evidence': ["Models, Paleo-observation, Contemporary observations"],
               'Key Attributes - Confidence: Existence of RS': ["Well established – Wide agreement in the literature that the RS exists"],
               'Key Attributes - Confidence: Mechanism underlying RS': ["Contested – Multiple proposed mechanisms, reasonable evidence both for and against different mechanisms"],
               'Links to other regime shifts': ["Impoundments, Shrub encroachment, Monsoon circulation, Forest to Savannas"]
              }
df_dict = pd.DataFrame(desert_dict)
RSD = RSD.merge(df_dict,how='outer')

## Converting '...Typical time scale' => 'CollapseTime' and '...Reversibility' => 'RecoveryTime'. This is the methodology of Barfuss et al. 2024
rev_nums = []
for x in RSD[['Key Attributes - Reversibility','Key Attributes - Typical time scale']].values:
    rev_nums.append(revstring_to_number(x[0],x[1]))
RSD['CollapseTime'] = RSD['Key Attributes - Typical time scale'].apply(tstring_to_number)
RSD['RecoveryTime'] = pd.Series(rev_nums)
RSD[['Regime Shift Name','RecoveryTime','CollapseTime','Key Attributes - Typical spatial scale']].head()

## Turning all regime shift names into lowercase for easier semantic cohesion
RSD['Regime Shift Name']=RSD['Regime Shift Name'].str.lower()

## Semantic Cohesion, i.e. making sure the 'Regime Shift' values from CLD align exactly with the 'Regime Shift Name' from RSD
# Manually making sure (similar) names are aligned
RSD=RSD.replace(to_replace='coniferous to deciduous boreal forest',value='coniferous to deciduous forest')
RSD=RSD.replace(to_replace='forest to savannas',value='forest to savanna')
df_edgelist=df_edgelist.replace(to_replace='kelps transitions',value='kelp transitions')
df_edgelist=df_edgelist.replace(to_replace='salt marshes to tidal flats',value='salt marsh to tidal flat')
RSD=RSD.replace(to_replace='marine food webs: community change and trophic level decline',value='marine foodwebs')
RSD=RSD.replace(to_replace='tundra to boreal forest',value='tundra to forest')
RSD=RSD.replace(to_replace='west antarctic ice sheet collapse',value='wais')
RSD=RSD.replace(to_replace='coastal marine eutrophication',value='marine eutrophication')
df_edgelist=df_edgelist.replace(to_replace='primary production arctic ocean',value='primary production in the arctic ocean')
df_edgelist=df_edgelist.replace(to_replace='floating plants',value='submerged to floating plants') # NB this will slightly change the order of regime shifts compared to Rocha et al 2018

RSD=RSD.replace(to_replace='indian summer monsoon',value='monsoon')
df_edgelist=df_edgelist.replace(to_replace='moonson',value='monsoon') # mispellings in Rocha et al 2018!

RSD=RSD.replace(to_replace='river channel position',value='river channel position change')
df_edgelist=df_edgelist.replace(to_replace='river channel change',value='river channel position change') # mispellings in Rocha et al 2018!

RSD=RSD.replace(to_replace='thermokarst lake to terrestrial ecosystem ',value='thermokarst lake to terrestrial ecosystem')
df_edgelist=df_edgelist.replace(to_replace='thermokarst lakes',value='thermokarst lake to terrestrial ecosystem')

## If you need to check what names are present in both uncomment the next two lines and print out the results. 
    #Can then do a symmetric difference (set1 ^ set2) to see what isn't in the intersection, 
    #or a set difference (set1 -set2) to see what set1 contains that isnt in set2. 
    #NB RSD has 2 regime shifts which are not in CLD
# rsd_set = set(sorted(RSD['Regime Shift Name'].unique()))
# cld_set = set(sorted(df_edgelist['Regime Shift'].unique()))

## Processing CLD edgelist into various NetworkX DiGraphs (directed graphs)
# turning pandas dataframe (edgelist) into a full causal loop diagram G_cld
# Nodes here are categorical variables, edges showing causal relationship
R = df_edgelist['Regime Shift'].unique() # producing an alphabetical ordering to the regime shifts in CLD

# G_cld is the ENTIRE causal loop diagram of all the 30 regime shifts 
G_cld = nx.from_pandas_edgelist(df_edgelist,source='Tail',target='Head',create_using=nx.DiGraph(),edge_attr=True)

# G_r_list is a dictionary with each value being the causal loop diagram of an individual regime shift 
    # e.g. G_r_list['wais'] is the cld for the wais regime shift
G_r_list = {r_name:nx.from_pandas_edgelist(df_edgelist[df_edgelist['Regime Shift']==r_name],
                                           source='Tail',target='Head',create_using=nx.DiGraph(),edge_attr=True) 
            for r_name in R}
V_r_list = {r:G_r_list[r].nodes() for r in R} # the corresponding dictionary of nodes

# G_r1r2_list is a dictionary of pairwise joint regime shift clds; 
    # e.g. G_r1r2_list['wais','desertification'] is the joint causal loop diagram between WAIS and Desertification
G_r1r2_list = {(r1,r2):nx.from_pandas_edgelist(df_edgelist[(df_edgelist['Regime Shift']==r1) | (df_edgelist['Regime Shift']==r2)],
                                               source='Tail',target='Head',create_using=nx.DiGraph(),edge_attr=True) 
            for r1 in R for r2 in R}

### Using the pairwise joint clds and the large cld to generate the three response variable matrices 
    # (i.e. three new adjacency matrices / graphs) as defined in Rocha et al 2018:
        # 1. Driver sharing (DS) - two regime shifts share an undirected edge if the both share drivers (nodes not in a cycle in individual clds)
        # 2. Domino effect (DE) - edge (r1,r2) exists if there exists a node n where n is a driver in the cld(r1) but not a driver in the cld(r2)
        # 3. Hidden feedbacks (HF) - number of new cycles generated by combining two individual clds
    # As well as a fourth matrix to use in NEcoPG:
        # 4. NEcoPG graph (NEco) - the weight of edge (r1,r2) is the number of simple directed paths FROM r1 TO r2

## N.B. "Drivers" are defined strictly in the code/methodology of Rocha et al 2018 as nodes that do not belong in any cycles. 
    ## The supplementary occassionally mentions drivers being able to impact feedbacks but not the other way round (implying there is a path from
    ## driver to at least one cycle, but not the other way round (i.e. not cycle -> driver). 
    ## IN OTHER WORDS, MISTAKE IN THE PAPER VS CODE!!
DS_list = [] 
DS_att = {} 
DE_list = []
DE_att = {}
HF_list = []
HF_att = {}
CP_list = []
CP_att = {}
for u in R:
    for v in R:
        if v != u:
            G_u = G_r_list[u]
            G_v = G_r_list[v]
            G_uv = G_r1r2_list[(u,v)]
            G_uv_intersection = nx.intersection(G_u,G_v)
            u_cycles = list(nx.simple_cycles(G_u))
            v_cycles = list(nx.simple_cycles(G_v))
            u_cycle_nodes = [n for cycle in u_cycles for n in cycle]
            v_cycle_nodes = [n for cycle in v_cycles for n in cycle]
            
            edges_uv = [] # NEcoPG edges generation
            for n1 in V_r_list[u]:
                for n2 in V_r_list[v]:
                    if n1 != n2 and nx.has_path(G_uv,n1,n2):
                        edges_uv.append((n1,n2))           
            if len(edges_uv) > 0:
                CP_list.append((u,v))
                CP_att[(u,v)] = len(edges_uv)

            ## Graphs as defined by Rocha et al 2018
            # Hidden feedbacks response variable edges generation
            hf = len(list(nx.simple_cycles(G_uv))) - len(u_cycles) - len(v_cycles)
            if hf > 0:
                HF_list.append((u,v))
                HF_att[(u,v)] = hf

            ds = []
            de = []
            for n in G_uv_intersection.nodes():
                if n not in v_cycle_nodes: # checking that n is a driver in G_v
                    if n not in u_cycle_nodes:
                        ds.append(n) # if n is a driver in G_u -> driver sharing
                    else:
                        de.append(n) # if n is not a driver in G_u -> domino effect
            if len(de) > 0:
                DE_list.append((u,v))
                DE_att[(u,v)] = len(de)
                
            if len(ds) > 0:
                DS_list.append((u,v))
                DS_att[(u,v)] = len(ds)

## From the four edgelists (and their respective edge_attributes dictionaries) produce the four new NetworkX graphs
G_CP = nx.DiGraph()
G_CP.add_nodes_from(R)
G_CP.add_edges_from(CP_list)
nx.set_edge_attributes(G_CP,CP_att,name='weight')
nx.set_edge_attributes(G_CP,{n:CP_att[n]/(len(V_r_list[n[0]])*len(V_r_list[n[1]])) 
                               for n in CP_att},
                       name='normalized weight')
A_CP = nx.to_numpy_array(G_CP,nodelist=sorted(R),weight='normalized weight')

G_DS = nx.Graph()
G_DS.add_nodes_from(R)
G_DS.add_edges_from(DS_list)
nx.set_edge_attributes(G_DS,DS_att,name='weight')
A_DS = nx.to_numpy_array(G_DS,nodelist=sorted(R),weight='weight')

G_DE = nx.DiGraph()
G_DE.add_nodes_from(R)
G_DE.add_edges_from(DE_list)
nx.set_edge_attributes(G_DE,DE_att,name='weight')
A_DE = nx.to_numpy_array(G_DE,nodelist=sorted(R),weight='weight')

G_HF = nx.DiGraph()
G_HF.add_nodes_from(R)
G_HF.add_edges_from(HF_list)
nx.set_edge_attributes(G_HF,HF_att,name='weight')
A_HF = nx.to_numpy_array(G_HF,nodelist=sorted(R),weight='weight')

R=sorted(R)

# Adding recovery and collapse times as node attributes
RSD = RSD[~RSD['Regime Shift Name'].isin(['common pool resource harvesting', 'invasive floating to invasive submerged plant dominance'])]
RSD = RSD.sort_values(by=['Regime Shift Name'])
RSD_node_atts = RSD[['Regime Shift Name','RecoveryTime','CollapseTime']].set_index('Regime Shift Name').to_dict('index')
RSD_indexmap = {'arctic benthos borealisation':'ABB', 
                'arctic sea-ice loss': 'ASI', 
                'bivalves collapse':'BIVS', 
                'bush encroachment':'BUSH',
                'coniferous to deciduous forest':'CODE',
                'coral transitions': 'REEF', 
                'desertification':'DSRT',
                'fisheries collapse': 'FISH',
                'forest to savanna': 'FOSA',
                'freshwater eutrophication': 'FETR',
                'greenland ice sheet collapse': 'GIS',
                'hypoxia': 'HYPX',
                'kelp transitions':'KELP',
                'mangroves transitions': 'MNGR',
                'marine eutrophication': 'METR',
                'marine foodwebs': 'MFDW',
                'monsoon': 'MONS',
                'peatland transitions': 'PEAT',
                'primary production in the arctic ocean': 'PPAC',
                'river channel position change': 'RVRC',
                'salt marsh to tidal flat': 'SMTF',
                'seagrass transitions': 'SGRS',
                'soil salinization': 'SOIL',
                'sprawling vs compact city': 'CITY',
                'steppe to tundra': 'STPT',
                'submerged to floating plants': 'FLOT',
                'thermohaline circulation': 'THLC',
                'thermokarst lake to terrestrial ecosystem': 'THKL',
                'tundra to forest': 'TUFO',
                'wais': 'WAIS'
               }
RSD['index'] = RSD['Regime Shift Name'].map(RSD_indexmap)
RSD = RSD.set_index('index')

nx.set_node_attributes(G_CP,RSD_node_atts)
nx.set_node_attributes(G_DS,RSD_node_atts)
nx.set_node_attributes(G_DE,RSD_node_atts)
nx.set_node_attributes(G_HF,RSD_node_atts)