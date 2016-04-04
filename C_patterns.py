from lingpy import *
import networkx as nx
import json
from sys import argv
from L_edit import pcdist, simple, editdist

def prepare_data():

    wl = Wordlist('D_wang-2006.tsv')
    restriction = [x[0] for x in csv2list('D_reconstructible-words.tsv')]

    wl.add_entries('characters', 'partial',lambda x: x)

    return wl, restriction

# define the three different models employed in the paper
model = 'simple'
if 'pcd' in argv:
    model = 'pcd'
elif 'edit' in argv:
    model = 'edit'
elif 'lump' in argv:
    model = 'lump'

# check for multiple or binary state
if 'multistate' in argv:
    state = 'multistate'
else:
    state = 'binary'

# load and create the dictionary
try:
    D = json.loads(open('D_data.json').read())
except:
    D = {}

# start data-preparation
wl, restriction = prepare_data()

# define the data that shall be produced
patterns = []
characters = []
matrices = []
concepts = []
factors = {}

# start the creation of the data
for concept in wl.concepts:
    if concept in restriction:

        print('? '+concept)

        # get the indices
        idxs = [k for k in wl.get_list(concept=concept, flat=True) if wl[k,'characters']]
        
        # get chars and taxa
        chars = [wl[k,'characters'] for k in idxs]
        taxa = [wl[k,'taxa'] for k in idxs]

        # get connected components
        states = sorted(set([wl[k,'characters'] for k in idxs]))

        if state == 'binary':
            if model == 'simple':
                comps = [[s] for s in states]
            elif model == 'lump':
                _G = nx.Graph()
                for i,s1 in enumerate(states):
                    for j,s2 in enumerate(states):
                        _G.add_node(s1)
                        if i < j:
                            if [x for x in s1 if x in s2]:
                                _G.add_edge(s1, s2)

                comps = list(nx.connected_components(_G))
                new_states = {}
                for c in comps:
                    this_char = sorted(c, key=lambda x: len(x),
                            reverse=False)[0]
                    for char in c:
                        new_states[char] = this_char
                states = sorted(set(new_states.values()))
                for i,c in enumerate(chars):
                    chars[i] = new_states[c]

                comps = [[s] for s in states]

        else:
            comps = [states]

        for comp in comps: 
            # get the taxon for a comp
            tmp = {}
            for char in comp:

                current_taxa = []
                for i in range(len(chars)):
                    if chars[i] == char:
                        try:
                            tmp[taxa[i]] += [char]
                        except KeyError:
                            tmp[taxa[i]] = [char]


            # now fill the rest with spaces
            pattern = []
            for taxon in wl.taxa:
                if taxon in tmp:
                    pattern += [tmp[taxon]]
                else:
                    pattern += [['-']]

            goon = True
            if pattern.count(['-']) >= len(wl.taxa)-1:
                goon = False

            if goon:
                print('-> '+concept)
                patterns += [pattern]
                concepts += [concept]
                characters += [sorted(comp)+['-']]

                freqs = {}
                for c in characters[-1]:
                    freqs[c] = 0
                    for p in pattern:
                        if c in p:
                            freqs[c] += 1

                # create the matrix
                matrix = [[0 for t in characters[-1]] for u in
                        characters[-1]]
                for i,sA in enumerate(characters[-1]):
                    for j,sB in enumerate(characters[-1]):
                        if i != j:
                            if 'pcd' in argv:
                                matrix[i][j] = pcdist(sA, sB)
                            elif 'edit' in argv:
                                matrix[i][j] = editdist(sA, sB)
                            elif 'simple' in argv:
                                matrix[i][j] = simple(sA, sB)
                            elif 'lump' in argv:
                                matrix[i][j] = 0 if [x for x in sA if x in sB] else 1
                            
                try:
                    factors[concept] += 1
                except KeyError:
                    factors[concept] = 1

                matrices += [matrix]

# correct for concepts 
for idx,(concept,matrix) in enumerate(zip(concepts,matrices)):
    factor = 1 / factors[concept]
    print(idx+1,concept, factor)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i < j:
                matrices[idx][i][j] = int(matrix[i][j] * factor + 0.5)
                matrices[idx][j][i] = int(matrix[j][i] * factor + 0.5)

D[model+'.'+state+'.matrices'] = matrices
D[model+'.'+state+'.chars'] = characters
D[model+'.'+state+'.patterns'] = patterns
D[model+'.'+state+'.concepts'] = concepts
D[model+'.'+state+'.taxa'] = wl.taxa

with open('D_data.json', 'w') as f:
    f.write(json.dumps(D, indent=2))

print('\n'.join(['* '+x for x in restriction if x not in concepts]))

print(model, state)

