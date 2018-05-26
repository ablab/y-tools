#! /usr/bin/env python

import argparse
import subprocess
from os import listdir, rmdir, mkdir
from shutil import rmtree
from math import log 
import os
import pandas as pd


import pygraphviz as gvz
from bokeh.io import show, output_file, save
from bokeh.plotting import figure
from bokeh.models import GraphRenderer, StaticLayoutProvider, Oval, HoverTool, \
    TapTool, BoxSelectTool, MultiLine, WheelZoomTool, PanTool, Plot, Range1d, Circle
import networkx as nx
from bokeh.models.graphs import NodesAndLinkedEdges, EdgesAndLinkedNodes

from holoviews.operation.datashader import datashade, bundle_graph
import holoviews as hv

def shms_pair_to_string(shms):
    V_shms, J_shms = shms
    if not isinstance(V_shms, str):
        V_shms = ''
    if not isinstance(J_shms, str):
        J_shms = '' 
    res = '\n'.join(["V_shms", '\n'.join(["\t" + x for x  in V_shms.split(";")]),
                     "J_shms", '\n'.join(["\t" + x for x  in J_shms.split(";")])]) 
    return res
def shms_to_string(shms):
    if not isinstance(shms, str):
        return ''
    res = '  |  '.join([x for x  in shms.split(";")]) 
    return res

def abundance_to_size(n):
    coeff = 1
    if n < 10: 
        coeff = 0.8
    elif n < 100:
        coeff = 1.6
    elif n < 1000:
        coeff = 2.4
    else:
        coeff = 3.2

    coeff = log(n+3, 4)
    return coeff*1.2, 1.2 #coeff*0.8

def draw_tree(antevolo_res_dir, tree_name, output_dir):
    tree_file = os.path.join(antevolo_res_dir, 'clonal_trees/', tree_name)
    print("drawing "+tree_file)
    vertices_file = os.path.join(antevolo_res_dir,'clonal_trees_vertices/', tree_name)
    vertex_to_depths = {}
    depth_to_vertices = {}
    vertex_to_mutations = {}
    edge_to_mutations = {}
    edges = []
    max_depth = 0
    vertices = {}
    clones = {}
    passed = set()

    vertices_df = pd.read_csv(vertices_file, sep="\t")
    tree_df = pd.read_csv(tree_file, sep="\t")

    for i, row in vertices_df.iterrows():
        clone_num = row["Clone_id"]
        clone_name = row["Clone_name"]
        clone_AA_seq = row["AA_seq"]
        clone_productive = row["Productive"]
        clone_left_anchor_AA = row["Left_CDR3_anchor_AA"]
        clone_right_anchor_AA = row["Right_CDR3_anchor_AA"]
        clones[clone_num] = [clone_productive, clone_AA_seq, clone_left_anchor_AA, clone_right_anchor_AA]
        clone_abundance = row["Size"]
        vertex_to_mutations[str(clone_num)] = (row["V_shms"], row["J_shms"])

        if clone_productive:
            clone_shape = 'ellipse'
        else:
            clone_shape = 'box'

            #clone_abundance = int(clone_name.split('_')[-1].split('|')[0])
        clone_width, clone_height = abundance_to_size(clone_abundance)
        if clone_name.split('_')[0] == 'fake':
            clone_color = 'magenta'
        else:
            clone_color = 'cyan'
        vertices[clone_num] = ''.join(['[label=',"\"" + str(clone_num)+'_'+clone_left_anchor_AA+clone_right_anchor_AA + "\"", 
                                       ', fixedsize=true, style=filled, fillcolor=', clone_color, 
                                       ', shape=', clone_shape, 
                                       ' width=', str(clone_width), ' height=', str(clone_height), ']'])

    for i, row in tree_df.iterrows():
        src_num = row["Src_id"]
        dst_num = row["Dst_id"]
        edge_type = row["Edge_type"]
        src_depth = row["Num_Src_SHMs"]
        dst_depth = row["Num_Dst_SHMs"]

        vertex_to_depths[src_num] = src_depth
        vertex_to_depths[dst_num] = dst_depth
        
        depth_to_vertices.setdefault(src_depth, set())
        depth_to_vertices[src_depth].add(src_num)
        depth_to_vertices.setdefault(dst_depth, set())
        depth_to_vertices[dst_depth].add(dst_num)
        
        edge_to_mutations[(str(src_num), str(dst_num))] = (row["V_shms"], row["J_shms"])

        edges.append([src_num, dst_num, edge_type, src_depth, dst_depth])   
        max_depth = max(max_depth, dst_depth)

        passed.add(src_num)
        passed.add(dst_num)

    fake_vertices = ['Depth_'+str(i) for i in range(max_depth+1)]
    DOT_OUTPUT_FILE_NAME = os.path.join(output_dir, tree_file.split('/')[-1]+".dot")
    with open(DOT_OUTPUT_FILE_NAME, 'w') as otp:
        otp.write("digraph "+'tree'+' {\n')
        otp.write("\tranksep=equally;")
        otp.write(''.join( ["\t{\n\t\tnode [shape=box]\n\t\t\n\t\t", ' -> '.join(fake_vertices), ";\n\t}\n\n"] ))
        for depth in depth_to_vertices:
            otp.write(''.join( ["\t{ rank = same;\n \t\t", 
                                "Depth_"+str(depth)+"; ", 
                                "; ".join(["\""+str(num)+"\"" for num in depth_to_vertices[depth]]),
                                ";\n\t};\n"] ))
        
        for v in vertices:
            if v in passed:
                otp.write(''.join(["\t","\""+str(v)+"\"", vertices[v],";\n"]))

        for edge in edges:
            src_num, dst_num, edge_type, src_depth, dst_depth = edge
            if edge_type == 'undirected' and src_depth == 0:
                continue

            if clones[src_num][1] == clones[dst_num][1]:
                edge_style = 'dotted'
            else:
                edge_style = 'filled'

            if clones[src_num][2] != clones[dst_num][2] or clones[src_num][3] != clones[dst_num][3]:
                edge_color = 'red'
            else:
                edge_color = 'black'

            dir_attr = ''
            if edge_type == 'reverse_directed':
                dir_attr = ', dir=both, arrowhead=none'
                src_num, dst_num = dst_num, src_num
            otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"",
                 " [color=", edge_color, ", style=", edge_style, dir_attr, "];\n"]))

        otp.write("}\n")
    subprocess.call(['dot', '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])
    return vertex_to_mutations, edge_to_mutations


def draw_bokeh_from_dot_file(input_filename, output_filename, vertex_to_mutations, edge_to_mutations):
    g = gvz.AGraph()

    g.read(input_filename)

    layout = g.layout()

    g.draw("{}.dot".format(output_filename), prog='dot', format="plain")

    g_with_layout = {}
    edges_layout_xs = []
    edges_layout_ys = []
    edge_starts = []
    edge_ends = []
    with open("{}.dot".format(output_filename)) as inp:
        for st in inp:
            arr = st.split()
            if arr[0] == 'node':
                g_with_layout[arr[1]] = (float(arr[2]), float(arr[3]))
            elif arr[0] == 'edge':
                edge_starts.append(arr[1])
                edge_ends.append(arr[2])
                edges_layout_xs.append([g_with_layout[edge_starts[-1]][0]] +\
                                       [float(coord) for i, coord in enumerate(arr[4:-2]) if i % 2 == 0] +\
                                       [g_with_layout[edge_ends[-1]][0]])
                edges_layout_ys.append([g_with_layout[edge_starts[-1]][1]] +\
                                       [float(coord) for i, coord in enumerate(arr[4:-2]) if i % 2 != 0] +\
                                       [g_with_layout[edge_ends[-1]][1]])



    max_x = max(v[0] for v in g_with_layout.values())
    max_y = max(v[1] for v in g_with_layout.values())
    plot = Plot(x_range=Range1d(-1.1,max_x + 1.1), y_range=Range1d(-1.1,max_y + 1.1), 
                 plot_width=1500, plot_height=900)
    plot.title.text = "Graph from tree"
    plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool(),
                   WheelZoomTool(), PanTool())
    graph = GraphRenderer()


    graph.node_renderer.data_source.add(list(g_with_layout.keys()), 'index')

    node_names = list(g_with_layout.keys())
    # node1_index = node_names.index(edge_starts[num])
    # node2_index = node_names.index(edge_ends[num])

    node_colors = ['blue']*len(g_with_layout.keys())
    # node_colors[node1_index] = 'red'
    # node_colors[node2_index] = 'black'
    graph.node_renderer.data_source.add(node_colors, 'color')
    graph.node_renderer.data_source.add(node_names, 'name')
    graph.node_renderer.glyph = Circle(radius=0.5, fill_color='color')
    graph.node_renderer.selection_glyph = Circle(radius=0.5, fill_color='black')
    graph.node_renderer.hover_glyph = Circle(radius=0.5, fill_color='red')

    node_v_shms = [shms_to_string(vertex_to_mutations[node][0])\
                            if node in vertex_to_mutations\
                            else '' for node in node_names]
    node_j_shms = [shms_to_string(vertex_to_mutations[node][1])\
                            if node in vertex_to_mutations\
                            else '' for node in node_names]
    graph.node_renderer.data_source.add(node_v_shms, 'V_shms')                    
    graph.node_renderer.data_source.add(node_j_shms, 'J_shms')                    

    graph.edge_renderer.data_source.data = dict(
        start=edge_starts,
        end=edge_ends)



    graph_layout = g_with_layout
    graph.layout_provider = StaticLayoutProvider(graph_layout=graph_layout)


    graph.edge_renderer.data_source.data['xs'] = edges_layout_xs
    graph.edge_renderer.data_source.data['ys'] = edges_layout_ys
    graph.edge_renderer.selection_glyph = MultiLine(line_color="blue", line_width=2)
    graph.edge_renderer.hover_glyph = MultiLine(line_color="orange", line_width=2)

    edge_v_shms = [shms_to_string(edge_to_mutations[(edge_start, edge_end)][0])\
                            if (edge_start, edge_end) in edge_to_mutations\
                            else '' for edge_start, edge_end in zip(edge_starts, edge_ends)]
    edge_j_shms = [shms_to_string(edge_to_mutations[(edge_start, edge_end)][1])\
                            if (edge_start, edge_end) in edge_to_mutations\
                            else '' for edge_start, edge_end in zip(edge_starts, edge_ends)]


    graph.edge_renderer.data_source.data['V_shms'] = edge_v_shms
    graph.edge_renderer.data_source.data['J_shms'] = edge_j_shms

    graph.selection_policy = NodesAndLinkedEdges()
    graph.inspection_policy = EdgesAndLinkedNodes()
    # bundled = bundle_graph(hv.Graph(graph))

    # plot.select_one(HoverTool).tooltips = [
    #     ('color', '@color'),
    #     ('name', '@name'),
    #     ('V_shms', '@V_shms'),
    #     ('J_shms', '@J_shms'),
    # ]
    plot.select_one(HoverTool).tooltips = [
        ('V_shms', '@V_shms'),
        ('J_shms', '@J_shms'),
    ]

    plot.renderers.append(graph)
    # plot.renderers.append(bundled)
    output_file("{}.html".format(output_filename))
    save(plot)
    # show(plot)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest = 'input', help="input dir with AntEvolo results", required=True)
    parser.add_argument("-o", "--output", dest = 'output', help="output dir", required=True)
    parser.add_argument("-s", "--strategy", dest = 'strategy',\
     help="'single' for specific tree (then -n TREE_FILE_NAME, 'topk' for a number of top-sized trees (then -k NUMBER_OF_TREES)", choices=['single', 'topk'] ,required=True)
    parser.add_argument("-k", "--trees_num", dest = 'k', type=int, help="number of top-size trees to draw")
    parser.add_argument("-n", "--name", dest = 'name', help="file name")
    args = parser.parse_args()

    if args.strategy == 'topk':
        trees = listdir(os.path.join(args.input, "clonal_trees/"))
        trees.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
        #for tree_file in trees[-int(args.k):]:
        #   print tree_file.split('_')[-1].split('.')[-2]
        
        try:
            listdir(args.output)
        except OSError:
            mkdir(args.output)
        for tree_name in trees[-int(args.k):]:
            vertex_to_mutations, edge_to_mutations = draw_tree(args.input, tree_name, args.output)
            tree_file = os.path.join(args.input, 'clonal_trees/', tree_name)
            DOT_OUTPUT_FILE_NAME = os.path.join(args.output, tree_file.split('/')[-1]+".dot")
            draw_bokeh_from_dot_file(DOT_OUTPUT_FILE_NAME, DOT_OUTPUT_FILE_NAME+".html", vertex_to_mutations, edge_to_mutations)
    else:
        draw_tree(args.input, args.name, args.output)

if __name__ == "__main__":
    main()