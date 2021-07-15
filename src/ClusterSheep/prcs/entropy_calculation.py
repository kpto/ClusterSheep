import logging
import math

from ClusterSheep.envr.session import get_session
from ClusterSheep.prcs.parallel.cluster_enrichment import enrich_clusters


try:
    session = get_session()
    clusters = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise


def calculate_entropy():
    _refresh_session()
    enrich_clusters()

    accumulated_clusters_entropy = 0
    num_idened_nodes = 0
    iden_stat = {}

    for cluster in clusters:
        accumulated_clusters_entropy, num_idened_nodes = accumulate_local_cluster_entropy(cluster, accumulated_clusters_entropy, iden_stat, num_idened_nodes)

    if len(iden_stat) == 0:
        logging.info("All spectra are unidentified, entropy cannot be calculated.")
        return None, None

    overall_cluster_entropy = accumulated_clusters_entropy / num_idened_nodes
    overall_peptide_entropy = 0

    for peptide, counts in iden_stat.items():
        peptides_entropy = 0
        total_count = sum(counts)
        for count in counts:
            probability = count / total_count
            peptides_entropy -= probability * math.log2(probability)
        overall_peptide_entropy += peptides_entropy * total_count / num_idened_nodes

    logging.info('Overall cluster entropy: {}'.format(overall_cluster_entropy))
    logging.info('Overall peptide entropy: {}'.format(overall_peptide_entropy))

    return overall_cluster_entropy, overall_peptide_entropy


def accumulate_local_cluster_entropy(cluster, accumulated_clusters_entropy, global_iden_stat, global_num_idened_nodes):
    iden_stat = cluster.identifications
    num_idened_nodes = sum(iden_stat.values())
    local_entropy = 0
    if num_idened_nodes > 0:
        counts_of_different_idens = iden_stat.values()
        for count in counts_of_different_idens:
            probability = count / num_idened_nodes
            local_entropy -= probability * math.log2(probability)

    accumulated_clusters_entropy += local_entropy * num_idened_nodes
    for iden, count in iden_stat.items():
        l = global_iden_stat.get(iden, None)
        if l is None:
            global_iden_stat[iden] = [count]
        else:
            l.append(count)
    return accumulated_clusters_entropy, global_num_idened_nodes + num_idened_nodes


def _refresh_session():
    global clusters
    clusters = session.clusters
    return


_refresh_session()


if __name__ == '__main__':
    pass