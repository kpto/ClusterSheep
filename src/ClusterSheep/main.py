# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 17:38:03 2017

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
from pathlib import Path
import sys

import ClusterSheep.envr.flags
import ClusterSheep.envr.session
import ClusterSheep.log
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
session = None
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def setup_environment():
    global session
    log_ = ClusterSheep.log.Log()
    log_.default()
    logging.info('Setting environment.')
    flags_ = ClusterSheep.envr.flags.Flags().parse_argv()
    ClusterSheep.envr.session.setup(flags_)
    session = ClusterSheep.envr.session.get_session()
    logging.info('Environment set.')
    log_.update()
    return


def main():
    try:
        if '--list-gpus' in sys.argv:
            from ClusterSheep.prcs.list_gpus import list_gpus
            list_gpus()
            exit()
        elif len(sys.argv) == 1:
            from ClusterSheep.prcs.help import print_help
            print_help()
            exit()

        setup_environment()

        if session.flags.print_session:
            print(session)
            exit()

        if session.flags.re_cluster:
            session.config.cg_finished.value = False
            session.config.cr_finished.value = False
            session.ranked_spectra = None
            session.clusters = None

        if session.flags.re_process:
            session.config.ii_finished.value = False
            session.config.id_finished.value = False
            session.config.rt_finished.value = False
            session.config.cg_finished.value = False
            session.config.cr_finished.value = False
            session.internal_index = None
            session.ranked_spectra = None
            session.clusters = None
            session.iden_lut = None

        from ClusterSheep.prcs.interactive_console import Checkpoint
        checkpoint = Checkpoint(session.flags.checkpoint, globals())

        if not session.config.ii_finished.value:
            checkpoint('Index Building')
            from ClusterSheep.prcs.parallel.index_building import build_index
            build_index()
        else:
            session.mount_internal_index()

        if not session.config.rt_finished.value:
            checkpoint('Rank Transformation')
            session.mount_internal_index()
            from ClusterSheep.prcs.parallel.rank_transformation import rank_transform
            rank_transform()
        else:
            if Path.cwd().joinpath(session.name + FILE_EXTENSION_RANKED_SPECTRA).exists():
                session.mount_ranked_spectra()

        if session.flags.rebuild_iden_lut or not session.config.id_finished.value:
            checkpoint('Identification Import')
            from ClusterSheep.prcs.parallel.identification_import import import_identification
            import_identification()
        else:
            if Path.cwd().joinpath(session.name + FILE_EXTENSION_IDEN_LUT).exists():
                session.mount_identification_lut()
            else:
                wrn_msg = 'Missing identification lookup table file.'
                logging.warning(wrn_msg)

        if not session.flags.preparation_only:

            if not session.config.cg_finished.value:
                checkpoint('Clustering')
                session.mount_internal_index()
                session.mount_ranked_spectra()
                from ClusterSheep.prcs.parallel.clustering import clustering
                clustering()
            if not session.config.cr_finished.value:
                checkpoint('Cluster Refinement')
                session.mount_clusters()
                from ClusterSheep.prcs.parallel.cluster_refinement import refine_cluster
                refine_cluster()
            if session.flags.no_saving:
                Path.cwd().joinpath(session.name + FILE_EXTENSION_RANKED_SPECTRA).unlink()
            session.mount_clusters()

            if session.flags.no_saving:
                checkpoint('Cluster Enrichment')
                session.mount_clusters()
                from ClusterSheep.prcs.parallel.cluster_enrichment import enrich_clusters
                enrich_clusters(True, session.config.gr_num_of_threads.value)
                checkpoint('Cluster Export')
                session.mount_clusters()
                from ClusterSheep.prcs.parallel.cluster_export import export_cluster
                export_cluster(Path.cwd().joinpath(session.name + '_exported_clusters'),
                               session.config.gr_num_of_threads.value)
                fmts = (FILE_EXTENSION_IDEN_LUT, FILE_EXTENSION_RANKED_SPECTRA, FILE_EXTENSION_CLUSTERS,
                        FILE_EXTENSION_RAW_CLUSTERS, FILE_EXTENSION_INTERNAL_INDEX, FILE_EXTENSION_SESSION)
                for f in list(Path.cwd().glob(session.name + '*')):
                    if f.is_file and f.suffix in fmts:
                        f.unlink()

        logging.info('All process finished.')
        if not session.flags.no_saving:
            session.save_session()

        if session.flags.stay_interactive:
            from ClusterSheep.prcs.cluster_viewing import cluster_viewer
            session.mount_internal_index()
            session.mount_clusters()
            cluster_viewer(globals())
            session.save_session()

        logging.info('Program ended.\n\n\n\n')
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            logging.info('Program terminated.\n\n\n\n')
        else:
            logging.exception('\nProgram ended unexpectedly. Logging traceback:\n'\
                              '==========TRACEBACK==========\n')
            print('\n\n\n\n')
            raise
# ====END OF CODE====
