# DrawParameterGraph.py
# Marcio Gameiro
# MIT LICENSE
# 2023-10-19

import DHN
import graphviz

def DrawFactorGraph(parameter_graph, node_index, node_color='lightblue',
                    ess_node_color='red', node_label='hex', node_size=None):
    """Draw logic factor graph of node given by node_index"""
    # Get list of hex codes
    hex_codes = parameter_graph.factorgraph(node_index)
    # Get list of vertices
    vertices = range(len(hex_codes))
    # Number of inputs
    n = len(parameter_graph.network().inputs(node_index))
    m = 1 # Number of thresholds
    ess_vertices = [v for v in vertices if DHN.essential(hex_codes[v], n, m)]
    # Check if two hex codes are adjacent
    adjacent = lambda u, v: DHN.isAdjacentHexcode(hex_codes[u], hex_codes[v])
    # Get list of edges
    edges = [(u, v) for u in vertices for v in vertices if adjacent(u, v)]
    # Create a vertex_name dictionary
    vertex_name = {}
    for v in vertices:
        vertex_name[v] = str(v)
    # Create a vertex_label dictionary
    vertex_label = {}
    for v in vertices:
        if node_label == 'hex':
            vertex_label[v] = hex_codes[v]
        elif node_label == 'index':
            vertex_label[v] = str(v)
        elif node_label == 'index:hex':
            vertex_label[v] = str(v) + ': ' + hex_codes[v]
        else: # Default ('hex')
            vertex_label[v] = hex_codes[v]
    # Node width parameters
    fixed_size = 'true'
    max_label_size = max([len(label) for label in vertex_label.values()])
    if node_size == None:
        node_size = max(0.4, max_label_size * 0.13)
        if max_label_size > 6:
            fixed_size = 'false'
    node_width = str(node_size)
    # Create a vertex_color dictionary
    vertex_color = {}
    for v in vertices:
        if v in ess_vertices:
            vertex_color[v] = ess_node_color
        else:
            vertex_color[v] = node_color
    # Now get a graphviz string for the graph
    graphviz_str = 'graph {' + '\n'.join(['"' + vertex_name[v] + '" [label="' + \
                    vertex_label[v] + '"; shape="circle"; fixedsize=' + fixed_size + \
                    '; width=' + node_width + '; style="filled"; fontsize=12; fillcolor="' + \
                    vertex_color[v] + '"];' for v in vertices]) + '\n' + \
                    '\n'.join(['"' + vertex_name[u] + '" -- "' + vertex_name[v]  + \
                    '";' for (u, v) in edges ]) + '\n' + '}\n'
    return graphviz.Source(graphviz_str)
