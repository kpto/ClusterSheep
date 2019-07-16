# -*- coding: utf-8 -*-
"""
Created on 16:59:37 27/11/2017

Author: Paul TO, Ka Po
Contact: kpto@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Division of Biomedical Engineering

Supervisor: Prof. Henry LAM, H. N.
Contact: kehlam@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Department of Chemical and Biomolecular Engineering

Desciption of this module:
"""

# ====BEGIN OF MODULE IMPORT====
import logging
import graph_tool.all as gt
import pickle
from pathlib import Path
from datetime import datetime
import os
from ast import literal_eval
import readline

from envr.session import get_session
from share.misc import generate_colors
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
cluster_viewer_banner = '''Input cluster id to view. For example:
viewer >>> 100

This command draws the cluster with id 100 in an interactive window.
Inside a interactive window, to select a spectrum, press "T" while pointing at it with the cursor.
Press "D" to plot the selected spectra. If there are two spectra selected, they will be plotted against each other.
You cannot select more than two spectra, when select, the second last is cleared and the last becomes the second last.
Press "I" to enable inverted identification. The identification of spectrum B is applied to spectrum A and vice versa.
Press "I" again to disable it.

To save image to disk, input "save {cluster_id} file={path} resolution={resolution} format={format}". For example:
viewer >>> save 100
    (Cluster is drawn and saved in working directory.
     The image file is named by date and time and the resolution is using default value 4000x4000)
viewer >>> save 100 file=graph_image.png
    (Cluster is drawn and saved in working directory with file name "graph_image.png")
viewer >>> save 100 file=/path/to/the/file
    (Cluster is drawn to the specified path. If the path is pointing to a directory,
     image file will be named by date and time)
viewer >>> save 100 resolution=1920x1080
    (Cluster is drawn and saved in working directory with resolution 1920x1080)
viewer >>> save 100 format=png
    (Cluster is drawn and saved as png image file, available formats are "ps", "pdf", "svg", and "png".
     If path specified includes an file extension, format argument will override the extension.)
viewer >>> save 100 file=~/graph_image format=png resolution=2000x2000

Input "enrich clusters" to append identification information to clusters.
This includes the major identification, the ratio of identified members to total and the number of identifications.
Enter "enrich clusters update" to update the identification information if the identification lookup table was modified.

Input "export {cluster_id} file={path}" to export a cluster to a text file. For example:
viewer >>> export 100
viewer >>> export 200 file=./export.txt
Source files and native ids of spectra in a cluster are listed line by line.
If file is not specified, the exported file will be named by date and time.

Input "quit" to quit.'''
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    internal_index = None
    iden_lut = None
    clusters = None
    iden_lut_cur = None
    clusters_cur = None
    iden_lut_conn = None
    clusters_conn = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class _Picked:

    def __init__(self, id_, entry, iden, type_):
        self.id_ = id_
        self.entry = entry
        self.iden = iden
        self.type_ = type_
        return

    def __str__(self):
        return str(self.id_)


class _Callback:

    def __init__(self):
        self.tagged_top = None
        self.tagged_bottom = None
        self.swap_iden = False
        self.verificative = True
        return

    def __call__(self, graphwidget, g, keyval, picked, pos, vprops, eprops):
        # key = T, tag vertex
        if keyval == 116:
            if 'spc' in g.vp:
                entry = g.vp['spc'][picked]
                entry.precursor_charge = 0
                type_ = 'spectrum'
            else:
                entry = session.internal_index[g.vp['iid'][picked]]
                type_ = 'index'
            identification = entry.get_identification()
            iden_string = identification.to_string() if identification else None
            self.tagged_bottom = self.tagged_top
            self.tagged_top = _Picked(picked, entry, identification, type_)
            logging.info('Picked: {}    Tagged pair: {}-{}'.format(picked, self.tagged_top, self.tagged_bottom))
            logging.info('Native id: {}    File: {}'.format(entry.native_id, entry.get_file_path()))
            logging.info('Precursor mass: {}    Charge: {}    Identification: {}'
                         .format(entry.precursor_mass, entry.precursor_charge, iden_string))
        # key = C, clear tags
        if keyval == 99:
            self.tagged_top = self.tagged_bottom = None
            logging.info('Tags cleared.')
        # key = D, plot spectrum
        if keyval == 100:
            spec_top = self.tagged_top.entry.get_spectrum() if self.tagged_top.type_ == 'index' else self.tagged_top.entry
            if self.tagged_bottom is None:
                spec_bottom = None
            else:
                spec_bottom = self.tagged_bottom.entry.get_spectrum() if self.tagged_bottom.type_ == 'index' else self.tagged_bottom.entry
                logging.info('Verificative ranked dot product: {}'.format(spec_top.verificative_ranked_dp(spec_bottom)))
                if self.swap_iden:
                    spec_top.override_iden = self.tagged_bottom.iden
                    spec_bottom.override_iden = self.tagged_top.iden
            spec_top.plot(against=spec_bottom, verificative=self.verificative)
        # key = I, swap identification
        if keyval == 105:
            self.swap_iden = False if self.swap_iden else True
            word = 'enabled' if self.swap_iden else 'disabled'
            logging.info('Identification swapping is {}.'.format(word))
        # key = V, trigger verificative plot
        if keyval == 118:
            self.verificative = False if self.verificative else True
            word = 'enabled' if self.verificative else 'disabled'
            logging.info('Verificative plotting is {}.'.format(word))
        return
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def cluster_viewer(globals_):
    _refresh_session()
    if iden_lut:
        iden_lut.connect()
    clusters.connect()
    _refresh_session()
    logging.info('Entering cluster viewer.')
    _read_history()
    print('\n' + '='*40)
    print(cluster_viewer_banner)
    if session.flags.dev_mode:
        print('\n' + '\033[1m\033[93mINPUT "python" TO ENTER PYTHON INTERACTIVE CONSOLE. DO ANYTHING YOU LIKE AT YOUR OWN RISK.\033[0m\033[0m')
    print('=' * 40 + '\n')
    while True:
        try:
            command = input('viewer >>> ')
        except EOFError:
            print('\n')
            break
        command = command.strip()
        if command == '':
            continue
        else:
            logging.debug('Viewer executed command: {}'.format(command))
        if command.strip() == 'quit':
            print('\n')
            break
        elif session.flags.dev_mode and command == 'python':
            from prcs.interactive_console import interactive_console
            _write_history()
            interactive_console(globals_)
            _read_history()
        elif command.startswith('enrich clusters'):
            _enrich_cluster(command)
        elif command.startswith('export'):
            _export(command)
        elif command.startswith('save'):
            _save(command)
        elif command.isdigit():
            command = int(command)
            if not _check_exists(command): continue
            _add_iden_stat(command)
            graph = _get_graph(command)[-1]
            draw_cluster_interactive(graph)
        else:
            logging.info('Invalid input.')
    logging.info('Leaving cluster viewer.')
    _write_history()
    return


def _read_history():
    history_file = Path.cwd().joinpath(session.name + FILE_EXTENSION_VIEWER_HISTORY)
    with history_file.open('w', encoding='utf-8') as fp:
        fp.write(session.viewer_hist)
    readline.read_history_file(str(history_file))
    return


def _write_history():
    history_file = Path.cwd().joinpath(session.name + FILE_EXTENSION_VIEWER_HISTORY)
    readline.write_history_file(str(history_file))
    with history_file.open(encoding='utf-8') as fp:
        session.viewer_hist = fp.read()
    history_file.unlink()
    return


def _save(string):
    args = string.split(' ')
    file = ''
    resolution = None
    format_ = None
    cluster_id = None

    for arg in args:
        if arg.isdigit():
            cluster_id = int(arg)
        elif arg.startswith('file='):
            file = arg.split('=')[1]
        elif arg.startswith('resolution='):
            resolution = arg.split('=')[1]
        elif arg.startswith('format='):
            format_ = arg.split('=')[1]
        elif arg == 'save':
            pass
        else:
            logging.info('Invalid input.')
            return

    if cluster_id is None:
        logging.info('Invalid input.')
        return

    valid_formats = ('.ps', '.pdf', '.svg', '.png')

    if file.startswith('~/'):
        file = Path.home().joinpath(file[2:]).absolute()
    elif file == '':
        file = Path.cwd()
    else:
        file = Path(file).absolute()
    if file.is_dir():
        file = file.joinpath('drawn_cluster_{}_'.format(cluster_id) + datetime.now().isoformat() + '.png')

    if not format_:
        format_ = file.suffix
        if format_ == '':
            format_ = '.png'
    else:
        format_ = '.' + format_

    if format_ not in valid_formats:
        logging.info('Invalid output format.')
        return
    file = file.parent.joinpath(file.stem + format_)

    if resolution:
        resolution = '({})'.format(resolution)
        resolution = resolution.replace('x', ',')
        resolution = resolution.replace('*', ',')
        try:
            resolution = literal_eval(resolution)
        except Exception:
            logging.info('Invalid resolution.')
            return
        if (type(resolution) is not tuple or len(resolution) != 2 or
           not all([type(x) == int for x in resolution]) or any([x <= 0 for x in resolution])):
            logging.info('Invalid resolution.')
            return
    else:
        resolution = (4000, 4000)

    if not _check_exists(cluster_id): return
    _add_iden_stat(cluster_id)
    graph = _get_graph(cluster_id)[-1]
    draw_cluster_save(graph, path=file, resolution=resolution)
    return


def _enrich_cluster(string):
    args = string.split(' ')
    num_of_threads = os.cpu_count()
    update = False

    for arg in args:
        if arg.isdigit():
            num_of_threads = int(arg)
        elif arg == 'update':
            update = True
        elif arg == 'enrich' or arg == 'clusters':
            pass
        else:
            logging.info('Invalid input.')
            return

    if num_of_threads < 1 or num_of_threads > os.cpu_count():
        logging.info('Invalid number of threads.')
        return

    from prcs.parallel.cluster_enrichment import enrich_clusters
    enrich_clusters(update, num_of_threads)
    _refresh_session()
    return


def _export(string):
    args = string.split(' ')
    cluster_id = None
    num_of_threads = os.cpu_count()
    all_ = False
    file = ''

    for arg in args:
        if arg.isdigit():
            cluster_id = int(arg)
        elif arg.startswith('file='):
            file = arg.split('=')[1]
        elif arg == 'all':
            all_ = True
        elif arg == 'export':
            pass
        else:
            logging.info('Invalid input.')
            return

    if not all_:
        if cluster_id is None:
            logging.info('Invalid input.')
            return
    else:
        num_of_threads = cluster_id if cluster_id else num_of_threads
        if num_of_threads < 1 or num_of_threads > os.cpu_count():
            logging.info('Invalid number of threads.')
            return

    if file.startswith('~/'):
        file = Path.home().joinpath(file[2:]).absolute()
    elif file == '':
        file = Path.cwd()
    else:
        file = Path(file).absolute()
    if file.is_dir():
        file = file.joinpath('exported_cluster_{}_'.format(cluster_id) + datetime.now().isoformat() + '.text')

    if not all_:
        if not _check_exists(cluster_id): return
        _add_iden_stat(cluster_id)
        _export_cluster(cluster_id, file)
    else:
        from prcs.parallel.cluster_export import export_cluster
        export_cluster(file, num_of_threads)
        _refresh_session()
    return


def _export_cluster(cluster_id, file):
    cluster = _get_graph(cluster_id)
    if cluster is None: return
    graph = cluster[-1]
    description = ('# cluster_id={}; num_nodes={}; num_edges={}; num_idens={}; '
                   'major_iden={}; iden_ratio={}; pre_mass_avg={}\n'.format(*(cluster[:-1])))
    files = session.ms_exp_files[session.internal_index.file_id[graph.vp['iid'].a]]
    native_ids = session.internal_index.native_id[graph.vp['iid'].a]
    if 'ide' in graph.gp:
        idens = graph.gp['ide']
        idens[-1] = 'None'
        iden_ids = graph.vp['ide'].a
        prb = graph.vp['prb'].a
        with file.open('a', encoding='utf-8') as fp:
            fp.write(description)
            for i in range(graph.num_vertices()):
                line = str(files[i]) + '\t' + str(native_ids[i]) + '\t' + idens[iden_ids[i]] + '\t' + str(prb[i]) + '\n'
                fp.write(line)
            fp.write('\n')
    else:
        with file.open('a', encoding='utf-8') as fp:
            fp.write(description)
            for i in range(graph.num_vertices()):
                line = str(files[i]) + '\t' + str(native_ids[i]) + '\n'
                fp.write(line)
            fp.write('\n')
    return


def _add_iden_stat(cluster_id, update=False):
    if not iden_lut:
        return
    graph = _get_graph(cluster_id)[-1]
    if 'ide' in graph.gp and not update:
        return
    num_vertices = graph.num_vertices()
    iid_arr = graph.vp['iid'].a
    ide = graph.new_vp('int', val=-1)
    prb = graph.new_vp('float', val=0.0)
    precursor_mass = np.empty(num_vertices, dtype=np.float32)
    ide_arr = ide.a
    prb_arr = prb.a
    idens = {}

    for i in range(num_vertices):
        entry = internal_index[int(iid_arr[i])]
        precursor_mass[i] = entry.precursor_mass
        identification = entry.get_identification()
        if identification and not identification.is_decoy:
            string = identification.to_tpp_string_integer() + '/' + str(identification.charge)
            if string not in idens:
                idens[string] = len(idens)
            ide_arr[i] = idens[string]
            prb_arr[i] = identification.probability
    pre_mass_avg = float(np.mean(precursor_mass))

    if len(idens) == 0:
        num_idens = 0
        iden_ratio = 0.0
        # clear original information in case of updating
        if 'ide' in graph.gp:
            graph.vp.pop('ide')
            graph.vp.pop('prb')
            graph.gp.pop('ide')
            graph.gp.pop('ord')
            pickled = pickle.dumps(graph)
            clusters_cur.execute('UPDATE "clusters" '
                                 'SET "num_idens"=?, "major_iden"=?, "iden_ratio"=?, "pre_mass_avg"=?, "pickled"=? '
                                 'WHERE "cluster_id"=?',
                                 (num_idens, None, iden_ratio, pre_mass_avg, pickled, cluster_id))
        else:
            clusters_cur.execute('UPDATE "clusters" '
                                 'SET "num_idens"=?, "iden_ratio"=?, "pre_mass_avg"=? '
                                 'WHERE "cluster_id"=?', (num_idens, iden_ratio, pre_mass_avg, cluster_id))
    else:
        num_idens = len(idens)
        graph.vp['ide'] = ide
        graph.vp['prb'] = prb
        idens = {v: k for k, v in idens.items()}
        graph.gp['ide'] = graph.new_gp('object')
        graph.gp['ide'] = idens
        iden_id, counts = np.unique(ide.a, return_counts=True)
        idx = np.where(iden_id == -1)[0]
        if len(idx) == 0:
            num_uniden = 0
        else:
            idx = idx[0]
            num_uniden = counts[idx]
            iden_id = np.append(iden_id[:idx], iden_id[idx+1:])
            counts = np.append(counts[:idx], counts[idx+1:])
        iden_id = iden_id[np.argsort(counts)[::-1]]
        graph.gp['ord'] = graph.new_gp('object')
        graph.gp['ord'] = list(iden_id)
        major_iden = idens[iden_id[0]]
        iden_ratio = (num_vertices - num_uniden) / num_vertices
        pickled = pickle.dumps(graph)

        clusters_cur.execute('UPDATE "clusters" '
                             'SET "num_idens"=?, "major_iden"=?, "iden_ratio"=?, "pre_mass_avg"=?, "pickled"=? '
                             'WHERE "cluster_id"=?',
                             (num_idens, major_iden, iden_ratio, pre_mass_avg, pickled, cluster_id))
    clusters_conn.commit()
    return


def _get_graph(cluster_id):
    query = clusters_cur.execute('SELECT * FROM "clusters" WHERE "cluster_id"=?',
                                 (cluster_id,)).fetchone()
    return query[:-1] + (pickle.loads(query[-1]),)


def get_graph(cluster_id=None):
    if cluster_id is None:
        logging.info('You need to input a cluster ID.')
        return
    query = clusters_cur.execute('SELECT * FROM "clusters" WHERE "cluster_id"=?',
                                 (cluster_id,)).fetchone()
    if not query:
        logging.info('Cluster with ID "{}" does not exist.'.format(cluster_id))
        return
    return query[:-1] + (pickle.loads(query[-1]),)


def _check_exists(cluster_id):
    is_exist = clusters_cur.execute('SELECT EXISTS(SELECT * FROM "clusters" WHERE "cluster_id"=?)',
                                    (cluster_id,)).fetchone()[0]
    if not is_exist:
        logging.info('Cluster with id "{}" does not exist.'.format(cluster_id))
        return False
    return True


def draw_cluster_interactive(graph, force_draw_edge=False):
    pos = gt.sfdp_layout(graph, p=1)
    if not force_draw_edge and graph.num_edges() > 10000:
        graph = gt.GraphView(graph, efilt=graph.new_ep('bool', val=False))
        logging.info('Number of edges exceeds 10000, edges are not drawn to ensure the stability.')
    fill_color = _fill_vertex_color(graph)
    gt.graph_draw(graph, pos=pos, vertex_size=9, vertex_fill_color=fill_color,
                  bg_color=(1, 1, 1, 1), key_press_callback=_Callback())
    return


def draw_cluster_save(graph, path=Path.cwd(), resolution=(4000, 4000)):
    pos = gt.sfdp_layout(graph, p=1)
    if path.is_dir(): path = path.joinpath(datetime.now().isoformat() + '.png')
    fill_color = _fill_vertex_color(graph)
    gt.graph_draw(graph, pos=pos, vertex_size=9, vertex_fill_color=fill_color,
                  bg_color=(1, 1, 1, 1), output_size=resolution, output=str(path))
    return


def _fill_vertex_color(graph):
    if 'ide' in graph.gp:
        ide = graph.vp['ide']
        order = graph.gp['ord']
        colors = generate_colors(len(order))
        color_map = {}
        for i in range(len(order)):
            color_map[order[i]] = colors[i]
        color_map[-1] = (0.0, 0.0, 0.0, 1.0)

        fill_color = graph.new_vp('vector<float>')
        for v in graph.vertices():
            fill_color[v] = color_map[ide[v]]
        return fill_color
    else:
        return 0, 0, 0, 1


def _refresh_session():
    global internal_index, iden_lut, clusters, iden_lut_cur, iden_lut_conn, clusters_cur, clusters_conn
    internal_index = session.internal_index
    iden_lut = session.iden_lut
    clusters = session.clusters
    if iden_lut:
        iden_lut_cur = iden_lut.cursor
        iden_lut_conn = iden_lut.connection
    if clusters:
        clusters_cur = clusters.cursor
        clusters_conn = clusters.connection
    return


_refresh_session()
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
