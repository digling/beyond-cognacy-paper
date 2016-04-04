from lingpy import *
from L_parsimony import *
from L_newick import *
import json
from sys import argv
import os

D = json.loads(open('D_data.json').read())

if not argv:
    print("Usage: python D_patterns.py [args]")
    
if 'summarize' in argv:

    data = csv2list('R_acr_all_results.txt', dtype=[str,float,float,int,str])

    txt = r"""  \tabular{|l||l|l|l||l|l|l||l|l|l|}
  \hline
  \multirow{2}*{\textbf{Model}} & 
  \multicolumn{3}{c||}{\textbf{\emph{Arbre}}} &
  \multicolumn{3}{c||}{\textbf{\emph{Shùxíngtú}}} &
  \multicolumn{3}{c|}{\textbf{\emph{Southern Chinese}}} \\\cline{2-10}
  & \textbf{Hits} & \textbf{Fails} & \textbf{Weight} & 
  \textbf{Hits} & \textbf{Fails} & \textbf{Weight} & \textbf{Hits} &
  \textbf{Fails} & \textbf{Weight}
  \\\hline\hline
  BINARY & {arbrebinary} & {sxtbinary} & {scbinary} \\\hline
  FITCH (multistate) & {arbrefitch} & {sxtfitch} & {scfitch}  \\\hline
  SANKOFF (multistate, weighted) & {arbresankoff} & {sxtsankoff} & {scsankoff} \\\hline
  DWST (multistate, weighed, directed) & {arbredwst} & {sxtdwst} & {scdwst} \\\hline
  \endtabular"""

    tmp = {
            'arbrepcd.multistate' : 'arbredwst',
            'arbreedit.multistate' : 'arbresankoff',
            'arbresimple.binary' : 'arbrebinary',
            'arbresimple.multistate' : 'arbrefitch',
            'shuxingtupcd.multistate' : 'sxtdwst',
            'shuxingtusimple.multistate' : 'sxtfitch',
            'shuxingtuedit.multistate' : 'sxtsankoff',
            'shuxingtusimple.binary' : 'sxtbinary',
            'southern_chinesepcd.multistate' : 'scdwst',
            'southern_chinesesimple.multistate' : 'scfitch',
            'southern_chineseedit.multistate' : 'scsankoff',
            'southern_chinesesimple.binary' : 'scbinary',
            }

    for line in data:
        key =  tmp[line[0]]
        goods = float(line[1])
        bads = float(line[2])
        alls = int(line[3])
        weight = line[4]
        hits = '{0:.2f} / {1:.2f}'.format(
                goods / alls,
                goods
                )
        fails = '{0:.2f} / {1:.2f}'.format(
                bads / alls,
                bads
                )

        rep_string = hits + '&' + fails + '&' + weight 
        txt = txt.replace('{'+key+'}',rep_string) #key, rep_string)

    with open('R_table.tex','w') as f:
        f.write(txt)
    
    import sys
    sys.exit()


# try to get a tree file
tree = False
for arg in argv:
    if arg.startswith('tree'):
        tree = LingPyTree(open(arg.split('=')[-1]).read())
        tree_name = arg.split('=')[-1][2:-4]

# search for matrix in argv
matrix = False
for arg in argv:
    if arg.startswith('matrix'):
        matrix = arg.split('=')[1]

# search for runs
runs = 15000
for arg in argv:
    if arg.startswith('runs'):
        runs = int(arg.split('=')[1])

# get the taxa from the matrix
taxa = D[matrix+'.taxa']

# this part of the script plots all data to html
if 'plot' in argv:
    
    folder = 'O_'+tree_name+'-'+matrix
    if not tree:
        raise ValueError("No tree specified!")

    if not os.path.isdir(folder):
        os.mkdir(folder)

    txt = '<h3> Parsimony Analysis for the {0} Model</h3>'.format(matrix) 
    txt += """<table class="table"><tr>
    <th>Number</th>
    <th>Weight</th>
    <th>Inferred Proto-Forms</th>
    <th>Concept</th>
    <th>Scenario</th>
    <th># (scenarios)</th></tr>
    """

    template = """
    <tr>
    <td>{0}</td>
    <td>{1}</td>
    <td>{2}</td>
    <td>{3}</td>
    <td>{4}</td>
    <td>{5}</td>
    </tr>"""
    
    total_weight = 0
    good = 0

    model_name = 'basic'
    
    # assign all labels for the data
    L = {}
    for node in tree.postorder:
        L[node] = {}
    
    # get the reconstructions
    all_reconstructions = []
    all_patterns = {}

    for idx,(p,m,c,cn) in enumerate(zip(D[matrix+'.patterns'], D[matrix+'.matrices'], 
        D[matrix+'.chars'], D[matrix+'.concepts'])):

        if cn not in ['knee','woman']:

            print('Analyzing {0} ({1})...'.format(idx+1, cn))
            debug = False
            if cn == 'woman':
                debug = True
            
            w,p,r = sankoff_parsimony(
                    p,
                    D[matrix+'.taxa'],
                    tree,
                    m,
                    c
                    )
            for pidx,pattern in enumerate(p):
                labels = {}
                _data = {}
                for t,c in pattern:
                    labels[tree[t]['label']] = '<b>['+c+']</b> <sup>'+tree[t]['label']+'</sup>'
                    _data[tree[t]['label']] = c
    
                tree.output('html',
                        filename=folder+'/pattern-'+model_name+'-'+str(idx)+'-'+str(pidx+1), labels=labels,
                        data=_data)
                try:
                    all_patterns[idx] += [pidx]
                except KeyError:
                    all_patterns[idx] = [pidx]

            
            # count the number of proposed reconstructions
            reconstructions = []
            for reconstruction in p:
                _tmp = dict(reconstruction)
                proto_form = _tmp[tree.root]
                reconstructions += [proto_form]
            # get the set
            pformset = sorted(set(reconstructions), key=lambda x:
                    reconstructions.count(x), reverse=True)
            pforms = ', '.join([x+' ({0:.2f})'.format(
                reconstructions.count(x) / len(reconstructions)
                ) for x in pformset])
            
            all_reconstructions += [reconstructions]

            txt += template.format(
                    idx,
                    w[0],
                    pforms,
                    cn,
                    ' '.join(['<a class="button" target="other" href="pattern-basic-{0}-{1}.html">{1}</a>'.format(
                            idx, pidx+1) for pidx in all_patterns[idx]]),
                    len(all_patterns[idx])
                    )
            total_weight += w[0]
    
    # create the style for the output-file
    style = """<style>
    .button {text-decoration: none; background: green; color: white; border: 1px solid gray;}
    .table {border: 2px solid black;} .table th {border: 2px solid gray;}
    .table td {border: 1px solid lightgray;}
    </style>"""

    # calculate innovations fore each node
    new_labels = {}
    new_data = {}
    new_labels['root'] = 'Old-Chinese'
    new_data['root'] = ''

    # table template
    ttemp = r"""<table><tr><th>Number</th><th>Proto-Form</th><th>Context</th><th>Parent Form</th><th>New Form</th></tr>"""
    
    # complete table
    txt += """<tr><td colspan="3">Total Weight</td><td>"""
    txt += str(total_weight)+'</td></td>'
    txt += '<td>'+'{0:.2f}'.format(0)+'</td>'
    txt += '<td>'
    txt += '<a target="other" href="pattern-'+model_name+'-summary.html" class="button">SHOW</a></td></tr></table>'

    with open(folder+'/index.html', 'w') as f:
        f.write('<html><head><meta charset="utf-8"</meta></head><body>'+style+txt+'</body></html>')

if 'acr' in argv:
    all_reconstructions = [] 
    all_weights = 0
    for idx,(p,m,c,cn) in enumerate(zip(D[matrix+'.patterns'], D[matrix+'.matrices'], 
        D[matrix+'.chars'], D[matrix+'.concepts'])):
        
        #print('Analyzing {0} ({1})...'.format(idx+1, cn))

        
        weight,chars = sankoff_parsimony_up(
                p,
                D[matrix+'.taxa'],
                tree,
                m,
                c,
                weight_and_chars = True
                )
        all_reconstructions += [chars]
        all_weights += weight

    och = dict(csv2list('D_old_chinese.csv'))
    #concepts = [x[1] for x in csv2list('data/words-wang.tsv') if x[2] == 'n']
    concepts = D[matrix+'.concepts']
    goods = 0
    bads = 0
    alls = 0

    if 'binary' in matrix:
        concept = ''
        tmp = []
        tmp2 = []
        for i,c in enumerate(concepts):
            if c != concept:
                concept = c
                tmp += [[x for x in all_reconstructions[i] if x != '-']]
                tmp2 += [c]
            else:
                tmp[-1] += [x for x in all_reconstructions[i] if x != '-']
        all_reconstructions = tmp
        concepts = tmp2

    with open('R_acr_'+matrix+'.'+tree_name+'.tsv', 'w') as f:
        f.write('NUMBER\tHIT\tCONCEPT\tPROTO\tSCORE\n')
        for idx,(c,r) in enumerate(zip(concepts, all_reconstructions)):

            # get the factor for the reconstructions
            if och[c] != '-':
                alls += 1

                factors = dict([(x, r.count(x) / len(r)) for x in set(r)])

                # get common items 
                commons = sorted(set([x for x in r if x in och[c].split(', ')]),
                        key=lambda x: factors[x], reverse=True)
                if not commons:
                    bads += 1
                    #print('! ','/'.join(sorted(set(r))), och[c])
                else:
                    good_score = sum([factors[x] for x in commons])
                    goods += good_score
                    bads += 1-good_score
                    #print('* ', commons[0], '{0:.2f}'.format(factors[commons[0]]), och[c])
                
                for com in commons:
                    f.write(str(idx+1)+'\t'+'+\t'+c+'\t'+com+'\t'+'{0:.2f}'.format(factors[com])+'\n')
                for com in [x for x in sorted(set(r)) if x not in commons]:
                    f.write(str(idx+1)+'\t'+'-\t'+c+'\t'+com+'\t'+'{0:.2f}'.format(factors[com])+'\n')


        print('Hits  (+):  ', '{0:.2f}'.format(goods / alls), '{0:.4f}'.format(goods))
        print('Fails (-): ', '{0:.2f}'.format(bads / alls),
                '{0:.4f}'.format(bads))
        print('Weight:', all_weights)
        f.write('\n---\n')
        f.write('+ {0:.2f}, {1:.2f}\n'.format(goods/alls,goods))
        f.write('- {0:.2f}, {1:.2f}\n'.format(bads/alls,bads))
        f.write('Weight: {0}'.format(all_weights))


        with open('R_acr_all_results.txt', 'a') as f:
            f.write(tree_name+matrix+'\t'+str(goods)+'\t'+str(bads)+'\t'+str(alls)+'\t'+str(all_weights)+'\n')



if 'innovations' in argv:

    transitions = {}
    for idx,(p,m,c,cn) in enumerate(zip(D[matrix+'.patterns'], D[matrix+'.matrices'], 
        D[matrix+'.chars'], D[matrix+'.concepts'])):

        print('Analyzing {0} ({1})...'.format(idx+1, cn))
        debug = False
        if cn == 'woman':
            debug = True
        
        w,p,r = sankoff_parsimony(
                p,
                D[matrix+'.taxa'],
                tree,
                m,
                c
                )

        
        tmp = dict(p[0])
        
        for node in tree.postorder[1:]:
            
            if not tree[node]['root']:
                previous_node = tree[node]['parent']
                next_char = tmp[node]
                previous_char = tmp[tree[node]['parent']]
                if next_char != previous_char:
                    try:
                        transitions[next_char] += [node]
                    except KeyError:
                        transitions[next_char] = [node]
    score = 0
    uniques = 0
    with open('R_transitions.'+matrix+'.'+tree_name+'.tsv', 'w') as f:
        for k,v in sorted(transitions.items(), key = lambda x: len(x[1])):
            for node in v:
                f.write(tree[node]['label'] + '\t'+node+'\t'+k+'\t'+str(len(v))+'\n')
            score += len(v)
            if len(v) == 1:
                uniques += 1
    print(score, uniques)
