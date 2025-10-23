#!/opt/miniconda3/envs/enan/bin/python

import json
import os
import numpy as np
from tabulate import tabulate
from collections import defaultdict

EH_TO_KCAL = 627.51
CAL_TO_J = 4.186
R = 8.314
T = 298.15

data = json.load(open(os.path.join(os.getcwd(), 'checkpoint.json')))



def boltz(array, T):
    return ((array-np.min(array))/(R*T)) / (np.sum((array-np.min(array))/(R*T)))

def get_table(data):
    rows = []

    for key, value in data.items():
        rows.append([key, *value.values()])

    return tabulate(rows, headers=['Protocol', 'Average', 'Number confs'], floatfmt='.10f')


protocols = defaultdict(list)
avs = defaultdict(dict)

for conf in data:
    for p in data[conf]['energies']:
        if data[conf]['energies'][p].get('G'): 
            protocols[p].append(data[conf]['energies'][p].get('G'))
        else:
            protocols[p].append(data[conf]['energies'][p].get('E'))

for key in protocols.keys(): 
    e = np.array(protocols[key]) * 1000 * CAL_TO_J
    pop = boltz(e, T)
    av = np.sum(e*pop)
    avs[key]['Average'] = av/(CAL_TO_J*1000*EH_TO_KCAL)
    avs[key]['Number confs'] = len(e)


print(get_table(avs))