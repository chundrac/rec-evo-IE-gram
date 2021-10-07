import itertools
from collections import defaultdict
import numpy as np
from numpy.random import uniform,normal
import re
from copy import copy
import os


tree = [l.rstrip().split('\t') for l in open('IE_nodes.csv','r')]

for i,l in enumerate(tree):
    tree[i] = l+['']*(max([len(m) for m in tree])-len(l))


branches = []

for i in range(len(tree[0])-1):
    k = []
    for j in range(len(tree)):
        if tree[j][i]!='':
            k.append(j)
    k.append(len(tree))
    for n in range(len(k)-1):
        sub = tree[k[n]:k[n+1]]
        for row in sub[1:]:
            if row[i+1]!='':
                branches.append((sub[0][i],row[i+1]))


daughters = defaultdict(list)

mothers = {}

for b in branches:
    mothers[b[1]] = b[0]


for b in branches:
    daughters[b[0]].append(b[1])


tips = [k for k in mothers.keys() if k not in daughters.keys()]


internal_nodes = [k for k in daughters.keys()]


nodes = tips+internal_nodes


dates = {}
for line in open('IE_dates.csv','r'):
#    dates[re.sub(r'\(|\)','',line.split('\t')[0])] = tuple([int(d) for d in line.split('\t')[1:3]])
    dates[line.split('\t')[0]] = [int(d) for d in line.split('\t')[1:3]]



for n in nodes:
    if n not in dates.keys():
        dates[n] = dates['_'.join(n.split('_')[:-1])]
        


prune_order = copy(tips)
nodes_to_order = copy(internal_nodes)
while len(nodes_to_order) > 0:
    for n in nodes_to_order:
        if len(daughters[n])==len([d for d in daughters[n] if d in prune_order]):
            prune_order.append(n)
            nodes_to_order.pop(nodes_to_order.index(n))


prune_order = prune_order[len(tips):]


branch_prune_order = []
daughters_prune_order = []
for n in prune_order:
    for d in daughters[n]:
        branch_prune_order.append([n,d])
        daughters_prune_order.append(d)


branch_num = []
for b in branch_prune_order:
    branch_num.append([nodes.index(b[0])+1,nodes.index(b[1])+1])


dates_ordered = []
for n in nodes:
    dates_ordered.append(dates[n])


for n in nodes:
    if n not in daughters.keys() or n.endswith('_advent'):
        dates[n][0] = dates[n][1]


np.random.seed(0)
branch_lengths = []
for i in range(100):
    dates_sampled = {}
    dates_sampled['Proto-Indo-European'] = uniform(dates['Proto-Indo-European'][0],dates['Proto-Indo-European'][1])
    for n in daughters_prune_order[::-1]:
#        dates_sampled[n] = normal(dates[n][0]+(dates[n][1]-dates[n][0])/2,25)
        dates_sampled[n] = uniform(max([dates[n][0],dates_sampled[mothers[n]]]),dates[n][1])
    lengths = []
    for b in branch_prune_order:
        if b[0] == b[1]+'_advent':
            lengths.append(1)
        else:
            lengths.append(abs(dates_sampled[b[1]]-dates_sampled[b[0]]))
    branch_lengths.append(lengths)



parent = [nodes.index(b[0])+1 for b in branch_prune_order]
child = [nodes.index(b[1])+1 for b in branch_prune_order]


# export to newick for visualization

def getdaughters(s):
#  if s in daughters.keys():
    return(daughters[s])


start = 'Proto-Indo-European'

def newick(s,i):
    daughterlist = getdaughters(s)
    if s not in mothers.keys():
        tree.insert(0,s+';')
    tree.insert(0,')')
    for d in daughterlist:
        tree.insert(0,d+':'+str(branch_lengths[i][daughters_prune_order.index(d)]))
        if daughterlist.index(d) != 0:
            tree.insert(1,',')
        if d in daughters.keys():
            newick(d,i)
    tree.insert(0,'(')


f = open('IE.newick','w')

for i in range(100):
    tree = []
    newick(start,i)
    print(''.join(tree),file=f)


f.close()