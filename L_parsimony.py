"""
Compute parsimony analyses.

Notes
-----

Code currently computes weighted parsimony. Required as input data are:

    * the patterns (the states in the leaves, passed as a list, multiple values
      are allowed and are interpreted as potential states of which the best
      states are then chosen)
    * taxa (the name of the languages, taxonomic units)
    * tree (the tree as a lingpy.cogent.object)
    * transitions (the matrix defining the transitions among the characters
    * characters: all characters that occur in the process

"""

import itertools
import random
import L_newick as nwk

def random_tree(taxa, branch_lengths=False):
    """
    Create a random tree from a list of taxa.

    Parameters
    ----------
    
    taxa : list
        The list containing the names of the taxa from which the tree will be
        created.
    branch_lengths : bool (default=False)
        When set to *True*, a random tree with random branch lengths will be
        created with the branch lengths being in order of the maximum number of
        the total number of internal branches.

    Returns
    -------
    tree_string : str
        A string representation of the random tree in Newick format.

    """
    # clone the list in order to avoid that lists used outside the function
    # suffer from modifications
    taxa_list = [t for t in taxa]

    random.shuffle(taxa_list)
    
    if not branch_lengths:
        while(len(taxa_list)  > 1):
            ulti_elem = str(taxa_list.pop())
            penulti_elem = str(taxa_list.pop())
            taxa_list.insert(0,"("+penulti_elem+","+ulti_elem+")")
            random.shuffle(taxa_list)
            
        taxa_list.append(";")
        return "".join(taxa_list)

    else:
        brlen_taxa_list = []
        nbr = 2*len(taxa_list)-3
        for taxon in taxa_list:
            brlen_taxa_list.append(str(taxon)+":"+'{0:.2f}'.format(random.uniform(1,nbr)))
        while(len(brlen_taxa_list) > 1):
            ulti_elem = str(brlen_taxa_list.pop())
            penulti_elem = str(brlen_taxa_list.pop())
            if len(brlen_taxa_list) > 0:
                brlen_taxa_list.insert(0,"("+penulti_elem+","+ulti_elem+")"+":"+'{0:.2f}'.format(random.uniform(0,nbr)))
            else:
                brlen_taxa_list.insert(0,"("+penulti_elem+","+ulti_elem+")")
            random.shuffle(brlen_taxa_list)
        brlen_taxa_list.append(";")
        return "".join(brlen_taxa_list)

def sankoff_parsimony_up(
        patterns, # the patterns in each taxonomic unit
        taxa, # the taxonomic units corresponding to the patterns
        tree, # the reference tree tree nodes in post-order
        transitions, # the transition matrix,
        characters, # the characters as they are provided in the transition matrix
        weight_only = False, # specify wether only weights should be returned
        weight_and_chars = False,
        debug = False
        ):
    """
    Carries out sankoff parsimony.

    Notes
    -----
    Think also of reducing parallel evolution by penalizing identical change
    patterns with extra weights! In this way, we car reduce the amount of
    parallel evolution, and also restrict the search space. Question is only
    how to handle this: Count all changes and penalize their occurrence?
    Problem is, we don't know in up-process, whether the change will really
    take place.
    """

    W = {}

    # get all characters

    # start iteration
    for node in tree.postorder:
        if debug:
            print('[D] Analyzing node {0} ({1} chars)...'.format(node,
                len(characters)))
        
        # name of node for convenience
        nname = node

        if tree[node]['leave']:

            W[nname] = {}
            for char in characters:
                if char in patterns[taxa.index(nname)]:
                    W[nname][char] = 0
                else:
                    W[nname][char] = 1000000
        else:
            
            W[nname] = {}

            # iterate over the states
            for nchar in characters:
                
                nscores = []
                
                # iterate over the children
                for child in tree[node]['children']:
                    cname = child

                    scores = []
                    for cchar in characters:
                        
                        # get the weight in child
                        wchild = W[cname][cchar]

                        # get the new weight due to transition process
                        wnew = wchild + transitions[
                                characters.index(nchar)
                                ][characters.index(cchar)]

                        # append to scores
                        scores += [wnew]

                    # get the minimal score for the char
                    smin = min(scores)

                    nscores += [smin]

                W[nname][nchar] = sum(nscores)
    
    if weight_only:
        return min(W[tree.root].values())
    
    if weight_and_chars:
        minw = min(W[tree.root].values())
        minchars = [x for x,y in W[tree.root].items() if y == minw]
        return minw,minchars
    
    return W
                    
def sankoff_parsimony_down(
        weights,
        patterns,
        taxa,
        tree,
        transitions,
        characters,
        debug = False
        ):
    
    # get the root
    root = tree.root

    # get the root chars
    smin = min(weights[root].values())

    # get the starting chars
    rchars = [a for a,b in weights[root].items() if b == smin]
    
    # prepare the queue
    queue = []
    for char in rchars:
        nodes = []
        for child in tree[tree.root]['children']:
            nodes += [child]
        queue += [([(nodes, tree.root, char)], [(tree.root, char)])]
    
    # prepare the scenarios which are written to output
    outs = []
    
    # start the loop
    while queue:

        nodes, scenario = queue.pop(0)

        if not nodes:
            outs += [scenario]
        else:
            
            if debug:
                print(len(queue))
            # get children and parent
            children, parent, pchar = nodes.pop()
            pidx = characters.index(pchar)

            # get the best scoring combination for scenario and children
            pscore = weights[parent][pchar]

            combs = itertools.product(*len(children) * [characters])
            
            for comb in combs:
                score = 0
                for i,char in enumerate(comb):
                    cidx = characters.index(char)
                    score += transitions[pidx][cidx]
                    score += weights[children[i]][char]
                
                if score == pscore:
                    new_nodes = [n for n in nodes]
                    new_scenario = [s for s in scenario]
                    
                    for child,char in zip(children,comb):
                        new_nodes += [(tree[child]['children'], child, char)]
                        new_scenario += [(child, char)]

                    queue += [(new_nodes, new_scenario)]
    return outs

def sankoff_parsimony(
        patterns,
        taxa,
        tree,
        transitions,
        characters,
        pprint=False,
        verbose = False,
        debug = False,
        ):
    if verbose:
        print('starting calculation')
    W = sankoff_parsimony_up(
            patterns,
            taxa,
            tree,
            transitions,
            characters,
            debug = debug
            )
    if verbose:
        print('starting backtrace')
    
    # get minimal weight
    smin = min(W[tree.root].values())
    weights = [b for a,b in W[tree.root].items() if b == smin]
    scenarios = sankoff_parsimony_down(
            W,
            patterns,
            taxa,
            tree,
            transitions,
            characters,
            debug
            )

    if pprint:
        tmp_tree = Tree(tree.newick)
        C = {}
        for k,v in tmp_tree.getNodesDict().items():
            C[str(v)[:-1]] = k

        for i,out in enumerate(scenarios):
            tr = tmp_tree.asciiArt()
            for k,v in out:
                target = v+len(C[k]) * '-'
                
                # get the nodes dict
                tr = tr.replace(C[k], target[:len(C[k])])
            print(tr)
            print(smin)
            print('')

    
    return weights, scenarios, W
