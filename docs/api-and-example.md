* [API documentation](#api-documentation)
    * [cluster_viewing](#api-documentation-cluster-viewing)
    * [cluster_enrichment](#api-documentation-cluster-enrichment)
    * [cluster_export](#api-documentation-cluster-export)
    * [Session](#api-documentation-session)
    * [Clusters](#api-documentation-clusters)
    * [Cluster](#api-documentation-cluster)
    * [Index](#api-documentation-index)
    * [Entry](#api-documentation-entry)
    * [Spectrum](#api-documentation-spectrum)
    * [IdentificationLUT](#api-documentation-identification-lut)
    * [RankedSpectra](#api-documentation-ranked-spectra)
* [Examples](#examples)


<a name="api-documentation"></a>
## API documentation

When entering Python interactive console, a set of modules and variables are automatically loaded for access of data. The loaded variables are as followed:

1. `cluster_viewing`, module [cluster_viewing](#api-documentation-cluster-viewing)
2. `cluster_enrichment`, module [cluster_enrichment](#api-documentation-cluster-enrichment)
3. `cluster_export`, module [cluster_export](#api-documentation-cluster-export)
4. `ses`, instance of [Session](#api-documentation-session)
5. `clu`, instance of [Clusters](#api-documentation-clusters)
6. `iid`, instance of [Index](#api-documentation-index)
7. `ide`, instance of [IdentificationLUT](#api-documentation-identification-lut)
8. `rks`, instance of [RankedSpectra](#api-documentation-ranked-spectra)

This documentation only includes methods and properties that are frequently used, some internal components are not documented.

Clusters are stored as `Graph` object from [graph_tool](https://graph-tool.skewed.de/) and it will be mentioned frequently, to better understand it or manipulate it in depth, see the documentation [here](https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.Graph).

<a name="api-documentation-cluster-viewing"></a>
### cluster_viewing

#### draw_cluster_interactive(`graph`, `force_draw_edge=False`)

This function is equivalent to entering a cluster ID in cluster viewer except that the input `graph` is a graph object instead of ID. It draws the cluster in an interactive window. If the number of edges of supplied cluster is more than 10000, edges will not be drawn to improve performance. To override this behaviour, set `force_draw_edge` to `True`.

#### draw_cluster_save(`graph`, `path=Path.cwd()`, `resolution=(4000, 4000)`)

This function is equivalent to command `save` in cluster viewer except that the input `graph` is a graph object instead of ID. It draws the cluster to a file specified by `path`. `resolution` is a tuple of horizontal and vertical pixels which specifies the resolution of the image.

<a name="api-documentation-cluster-enrichment"></a>
### cluster_enrichment

After clustering, clusters do not contain information of identification of spectra within. Cluster enrichment appends such information and internalize it in clusters. When a cluster is visualized, this process is first done so that nodes can be colored by identifications. Since this process involves random access of identification database, it takes a considerable amount of time when applied on many clusters.

#### enrich_clusters(`update=False`, `num_of_threads=os.cpu_count()`)

Enrich all clusters. Clusters that are already enriched are skipped if `update` is `False`. This method utilizes multi-processing to improve performance. The number of processes is specified by `num_of_threads` or all logical CPUs if it is not specified.

#### enrich_one(`cluster_id`, `update=False`)

Enrich one cluster. If the supplied cluster is already enriched and `update` is `False`, no calculation will be done. 

<a name="api-documentation-cluster-export"></a>
### cluster_export

#### export_clusters(`file`, `num_of_threads=os.cpu_count()`)

Export all clusters to a text file specified by `file` which expects a `Path` object. This method utilizes multi-processing to improve performance. The number of processes is specified by `num_of_threads` or all logical CPUs if it is not specified.

#### export_one(`file`, `cluster_id`)

Export one cluster to a text file specified by `file` which expects a `Path` object. If the file already exists, the cluster will be appended to it. 

<a name="api-documentation-session"></a>
### Session

#### name

Session name, specified by `--name=` argument.

#### creation_time

The time when this session is created.

#### config

Configuration used in this clustering session, specified by `--config=`. All parameters are printed when it is printed with `print()`.

#### ms_exp_files

Array of `Path` objects pointing to all MS experiment files.

#### iden_files

Array of `Path` objects pointing to identification files.

<a name="api-documentation-clusters"></a>
### Clusters

When iterated, it fetches row by row from the SQLite database and yields [Cluster](#api-documentation-cluster) objects. To prevent high RAM consumption, the `graph` property of yielded objects is lazy. The graph data is fetched from database when the property is accessed.

#### file_path

File path to the SQLite database storing all clusters.

#### get_cluster(`cluster_id`)

Returns [Cluster](#api-documentation-cluster) object.

#### exists(`cluster_id`)

Check does a cluster exist.

<a name="api-documentation-cluster"></a>
### Cluster

#### num_of_nodes

Number of nodes in this cluster.

#### num_of_edges

Number of edges in this cluster.

#### num_of_identifications

Number of identifications in this cluster.

#### major_identification

The identification shared by the most identified spectra.

#### identified_ratio

The ratio of identified spectra. 1 if all spectra have identification found in imported identification files.

#### average_precursor_mass

Average value of precursor mass of spectra in this cluster.

#### graph

`Graph` object from graph_tool. This property is lazy, the pickled object is loaded when it is accessed for the first time. The graph contains a set of property maps listed below:

1. Vertex property map `iid`, internal ID of spectra.
2. Edge property map `dps`, dot products of spectra pairs of edges.

If the graph is enriched, it contains additional property maps listed below:

3. Vertex property map `ide`, identifications of spectra, an integer being the key to access the string value stored in graph property map `ide` . -1 if a spectrum has no identification.
4. Vertex property map `prb`, probability of identification.
5. Graph property map `ide`, Python dictionary of all identifications keyed by integers.

For detail about property map, see [here](https://graph-tool.skewed.de/static/doc/quickstart.html#sec-property-maps).

<a name="api-documentation-index"></a>
### Index

`Index` is an array of lightweight representation of spectra. An [entry](#api-documentation-entry) of a spectrum can be fetched with `index_instance[internal_id]`.

<a name="api-documentation-entry"></a>
### Entry

`Entry` is a pointer of a spectrum.

#### get_file_path()

Get the path to the MS experiment file containing the spectrum. Return `Path` object.

#### get_identification()

Search identification database and get the identification of the spectrum, return `None` if it is not identified.

#### get_spectrum()

Return [Spectrum](#api-documentation-spectrum) object. This method reads the MS experiment file and get the peaks of the spectrum out of it.

<a name="api-documentation-spectrum"></a>
### Spectrum

`Spectrum` contains peak information of a spectrum.

#### mz

Array of MZ values of peaks.

#### intensity

Array of intensity values of peaks.

#### ranked

A `boolean` of is this spectrum rank-transformed.

#### clip(`mz_range=None`, `copy=True`)

Remove all peaks that are outside the specified MZ range. `mz_range` is a tuple of lower limit and upper limit. If it is `None`, the method uses the value in session configuration. Return a new `Spectrum` object if `copy` is `True`. Otherwise, the removal is applied on the current instance, and the current instance is returned.

#### remove_precursor(`removal_range=None`, `true_precursor_mass=None`, `copy=True`)

Remove all peaks near the precursor mass. `removal_range` is a tuple of lower addition and upper addition. For example, if (-20, 20) is used and the precursor mass is 400, peaks within 380-420 will be removed. `true_precursor_mass` is a boolean telling the method should the precursor mass be divided by the precursor charge. If `removal_range` or `true_precursor_mass` is `None`, the method uses the value in session configuration. Return a new `Spectrum` object if `copy` is `True`. Otherwise, the removal is applied on the current instance, and the current instance is returned.

#### rank_transform(`num_of_peaks=None`, `bins_per_th=None`, `copy=True`, `normalize=True`)

Rank-transform the spectrum. `num_of_peaks` specifies the number of the highest peaks extracted. `bins_per_th` specifies the resolution of bins by the number of bins per Thomson in MZ axis. For example, bin size will be 0.25Th if `bins_per_th` is 4. `normalize` decides should the ranks of peaks be normalized to become intensities. If `num_of_peaks` or `bins_per_th` is `None`, the method uses the value in session configuration. Return a new `Spectrum` object if `copy` is `True`. Otherwise, the removal is applied on the current instance, and the current instance is returned.

#### ranked_dot_product(`target`)

Calculate the dot product of two rank-transformed spectra. `target` expects another `Spectrum` object. If either one of spectra is not rank-tranformed, `None` is returned instead.

#### plot(`against=None`, `x_limit=(0, 2000)`, `verificative=True`, `show_iden=True`, `highlight=True`)

Plot the spectrum alone or against another spectrum specified by `against` which expects a `Spectrum` object. `x_limit` is a tuple of lower end and upper end of a range limiting the range of MZ axis being shown. `verificative` decides should the spectrum be processed according to the session configuration first, which includes clipping, removing precursor peaks and rank-transformation. If `show_iden` is `True`, identification database will be searched and if the spectrum has identification, the string representation of the identification will be printed in the plotted figure. If `highlight` is `True` and the spectrum has identification, the spectrum will be compared with a theoretical spectrum of the identification. Peaks matching ion peaks in the theoretical spectrum will be colored by the ion type.

<a name="api-documentation-identification-lut"></a>
### IdentificationLUT

#### file_path

File path to the SQLite database storing all imported identifications.

#### get_identification(`ms_exp_file_name`, `native_id`):

Get the identification by MS experiment file name specified by `ms_exp_file_name` and the scan number specified by `native_id`. `ms_exp_file_name` is the name of the file without parent path and the extension. For example, if the path of a file is `/path/to/file.mzXML`, the supplied name should be `file`. Return `None` if no identification is found.

<a name="api-documentation-ranked-spectra"></a>
### RankedSpectra<a name="api-documentation-s

All rank-transformed spectra that was used to compute dot products. This is the data uploaded to GPU.

#### mz

2D array. The first dimension is spectra and the second dimension is the MZ values of peaks of each spectrum.

#### intensity

2D array. The first dimension is spectra and the second dimension is the intensity of peaks of each spectrum.

#### num_of_peaks

An integer of the number of peaks in each spectrum. All spectra must have the same number of peaks.

