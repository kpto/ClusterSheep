![mascot-with-name](https://github.com/kpto/ClusterSheep/raw/release/docs/mascot-with-name.svg)

<div align="center">
    <div>
        <a href="https://dev.azure.com/kpto/ClusterSheep/_build?definitionId=4">
         <img src="https://img.shields.io/azure-devops/build/kpto/40d98952-bda4-49c3-b81b-d1f9debdfae7/4/release" /></a>
        <a href="https://dev.azure.com/kpto/ClusterSheep/_build?definitionId=3">
         <img src="https://img.shields.io/azure-devops/tests/kpto/40d98952-bda4-49c3-b81b-d1f9debdfae7/3" /></a>
        <a href="https://test.pypi.org/project/ClusterSheep">
         <img src="https://img.shields.io/pypi/v/ClusterSheep" /></a>
        <a href="https://hub.docker.com/r/kpto/clustersheep">
         <img src="https://img.shields.io/badge/docker-kpto%2Fclustersheep-099cec" /></a>
    </div>
</div>

---

* [Features](#features)
* [Overview](#overview)
* [Prerequisites](#prerequisites)
    * [Base](#prerequisites-base)
    * [Others](#prerequisites-others)
* [Installation](#installation)
    * [Docker](#installation-docker)
    * [PyPI](#installation-pypi)
* [Flow](#flow)
* [Usage](#usage)
    * [Quick start](#usage-quick-start)
    * [Cluster viewer](#usage-cluster-viewer)
        * [Visualization](#usage-cluster-viewer-cluster-visualization)
        * [Spectrum plotting](#usage-cluster-viewer-spectrum-plotting)
        * [Export](#usage-cluster-viewer-export)
    * [Advanced](#usage-advanced)
        * [Configuration](#usage-advanced-configuration)
        * [Material preparation](#usage-advanced-material-preparation)
        * [Python interactive console](#usage-advanced-python-interactive-console)
    * [Cluster refinement](#usage-cluster-refinement)
* [Limitation](#limitation)
* [Future](#future)
* [Bugs or requests](#bugs-or-requests)
* [Thesis](#thesis)


<a name="features"></a>
## Features

* Fast: CUDA accelerated, GPU (GTX 1070) performs pairwise similarity computation ~45 times faster than CPU (i7 6700K).
* Visualization: Powered by [graph-tool](https://graph-tool.skewed.de) and [matplotlib](https://matplotlib.org/), clusters are intuitively visualized using force directed drawing. Nodes can be picked to be plotted and peaks of peptide fragment ions are colored if the spectrum is identified. Detail on section [Cluster viewer](#usage-cluster-viewer).
* Multi-GPUs: Clustering with multiple GPUs is supported with almost no performance loss (>90% efficiency).
* Big data ready: ClusterSheep is written with big data in mind from the beginning, memory-map is heavily used to reduce memory usage and allow handling of large data set.
* Easy post-processing: Clusters and peptide identifications are stored in [SQLite](https://www.sqlite.org/index.html) database file for easy access. Output data can be manipulated easily with the built-in Python interactive console. Detail on section [Python interactive console](#usage-advanced-python-interactive-console).
* Traceability: Logging is important for research, ClusterSheep aggressively logs everything of a clustering session from its beginning to its last use, including configuration and all commands executed in cluster viewer and Python console.


<a name="overview"></a>
## Overview

ClusterSheep is a CUDA accelerated clustering software of MS2 spectra data generated from proteomics. ClusterSheep is inspired by [SpectraST](http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST), a software that builds spectral library by merging similar spectra together to form reference spectra. ClusterSheep works similarly to SpectraST, ClusterSheep groups spectra by computing similarities but with the following key differences.

1. ClusterSheep doesn't combine spectra, it outputs the resulted groups which are called clusters and also the pairwise similarity scores of similar spectra, represented as edges.

2. Unlike SpectraST and many existing clustering software which avoid pairwise similarity computation by producing a middle spectrum and compare candidate spectra with it only, ClusterSheep genuinely compute similarities of all pairs of spectra where the precursor mass difference of the pair is within the tolerance (1.1 by default).

3. ClusterSheep is designed to run on GPU instead of CPU. The computing power provided by GPUs makes pairwise computation possible within an acceptable time.

ClusterSheep is designed to improve the [spectral library building](https://en.wikipedia.org/wiki/Peptide_spectral_library) process. By genuinely computing pairwise similarities, the resulted cluster structure may reveal bridging (such as chimeric spectrum) of sub-clusters that have long been incorrectly merged during library building. Also, spectral clustering provides a middle step of updatable library building. New experiment data can be incorporated into existing clusters easily without re-computing old data. Though this has not yet been implemented in ClusterSheep.

ClusterSheep is written by Paul TO and supervised by [Prof. Henry LAM](https://facultyprofiles.ust.hk/profiles.php?profile=henry-hei-ning-lam-kehlam), the author of SpectraST.

<a name="prerequisites"></a>
## Prerequisites

<a name="prerequisites-base"></a>
### Base
* [Ubuntu](https://ubuntu.com) 18.04 or above *
* [Maxwell or above](https://en.wikipedia.org/wiki/List_of_Nvidia_graphics_processing_units#GeForce_900_series) CUDA capable graphics card  †
* [NVIDIA Driver](https://www.nvidia.com/Download/index.aspx)
* [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
* [graph-tool](https://graph-tool.skewed.de)
* [matplotlib](https://matplotlib.org)
* [pycuda](https://pypi.org/project/pycuda)
* [pyopenms](https://pypi.org/project/pyopenms)

<a name="prerequisites-others"></a>
### Others
* [Docker](https://docs.docker.com/engine/install/ubuntu) and [NVIDIA Container Toolkit](https://github.com/NVIDIA/nvidia-docker) (for installation via docker image) ‡
* [Cython](https://pypi.org/project/Cython) (for manually build from source code)

<br/>
<br/>

\* Any Linux distribution that can install [graph-tool](https://graph-tool.skewed.de) can run ClusterSheep, but installation on Ubuntu 18.04 is the easiest.

† Cards before Maxwell should run as well but they are not tested. Also, ClusterSheep can run without a graphics card, it fallbacks to CPU if no GPU is available.

‡ NVIDIA Container Toolkit works with Docker 19.03 or above. If the docker installed is below 19.03, [nvidia-docker2](https://github.com/nvidia/nvidia-docker/wiki/Installation-(version-2.0)) (deprecated) should be installed instead.

<a name="installation"></a>
## Installation

<a name="installation-docker"></a>
### Docker

This is the easiest way to install ClusterSheep, the docker image has already included everything ClusterSheep needs to run. You only need to install the NVIDIA Driver, Docker and NVIDIA Container Toolkit. Since the image is based on a [CUDA 9.2 image](https://hub.docker.com/layers/nvidia/cuda/9.2-devel-ubuntu18.04/images/sha256-42a0669c248d17225c3e336df2835cf017b299345c6431b3e364f37dae7645b2?context=explore), make sure that the version of installed NVIDIA Driver is >= 396.26. If you are installing on Ubuntu, the NVIDIA Driver can be installed via the built-in driver manager. With the above installed, the ClusterSheep image can be pulled with the following command:

```
docker pull kpto/clustersheep
```

<br/>
<br/>

It is recommended to start a container with following command:

```
docker run -ti --rm --gpus all -u user -w /home/user -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix kpto/clustersheep
```

`--gpus all` allows the container to access the GPU, `-u user` prevents running ClusterSheep as root, `-w /home/user` set the initial working directory to be an user space directory and `-e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix` enables X window rendering on host's X server. To access files of the host, you need one more `-v` option to map the desired directory. For example, `-v ~/Documents:/home/user/Documents` maps your Documents folder to the Documents folder inside the container. See [reference of docker run](https://docs.docker.com/engine/reference/run/#volume-shared-filesystems) for more information.

<a name="installation-pypi"></a>
### PyPI

ClusterSheep is available via PyPI, with all [basic dependencies](#prerequisites-base) installed, ClusterSheep can be installed with the following command:

```
pip install ClusterSheep
```

<a name="flow"></a>
## Flow

![flow](https://github.com/kpto/ClusterSheep/raw/release/docs/flow.svg)

The above is the flow of ClusterSheep. To begin with a new session, you need to name the session and supply MS experiment files and identification files (left side of the flow). ClusterSheep supports `mzXML`, `mzML` and `pep.xml` from either [PeptideProphet™](http://peptideprophet.sourceforge.net/) or InterProphet.

ClusterSheep will run clustering on provided data and produce a finished session, log and intermediate files. All files are named by the session name. A full clustering process contains following steps:

1. Create a new clustering session, store the session in a `cssess` file and create a logging file `cslogg`. All logging in the future will be written into it.
2. Build an index of all spectra which is stored in a `csindx` file.
3. Clarify (binning, precursor peaks removal and mz range clipping) and [rank-transform](https://pubmed.ncbi.nlm.nih.gov/24115759/) all spectra and store the transformed spectra in a `csrksp` file.
4. If any identification file is supplied, import identifications with probability equal or higher than the threshold and store them in a `csiden` file for future searching.
5. Pairwise similarity computation of all rank-transformed spectra where the precursor mass difference of the spectrum pair is within the tolerance. If the similarity score is higher than the threshold, an edge with the score are recorded. This step is accelerated by GPU.
6. Searching connected components from edges generated above. Each cluster is presented as a graph object of graph-tool. Graphs are stored in a `csclut` file.
7. Cluster refinement, optional, detail on section [Cluster refinement](#usage-cluster-refinement).

Both `csiden` and `csclut` are just [SQLite](https://www.sqlite.org/index.html) database files, they can be opened by any SQLite browser, for example, [DB Browser for SQLite](https://sqlitebrowser.org/).

Then, you can visualize the result by loading a finished session (right side of the flow). See section [Cluster viewer](#usage-cluster-viewer) for detail.

<a name="usage"></a>
## Usage

After installation, an entry point should be created so that ClusterSheep can be executed by simply executing command `clustersheep`. ClusterSheep assumes the current working directory to be the output directory, all files of a new clustering session will be generated on the directory where you start ClusterSheep.

Here are some example commands of auxiliary functions.

<br/>

```
clustersheep
```

With no parameter, ClusterSheep outputs a list of all available parameters. A configuration template file named `config_template` will be generated. This file contains all available parameters of ClusterSheep and the detail explanation of each setting.

<br/>

```
clustersheep --version
```

Which prints the version of the installed ClusterSheep.

<br/>

```
clustersheep --list-gpus
```

Which lists all available GPUs that can be used by ClusterSheep.

<a name="usage-quick-start"></a>
### Quick start

ClusterSheep accepts parameters with following formatting:

```
clustersheep --flag-type-option --value-type-option=value path-one path-two path-three ...
```

All options start with two hyphens `--`, with the value connected by an equal symbol `=` if the option expects a value. All non-option parameters in the end are seen as paths to MS or identification files or folders containing them. If the path is a folder, ClusterSheep searches supported files within it but it does not search sub-folders.

To cluster one MS experiment file, run ClusterSheep with `--name=` option as below:

```
clustersheep --name=mysession /path/to/the/ms-experiment-file.mzML
```

ClusterSheep will create a new clustering session named `mysession`. Files will be generated on the current working directory and are all named by the session name.

If you also have one identification file you want to import, execute below:

```
clustersheep --name=mysession /path/to/the/ms-experiment-file.mzML /path/to/the/identification-file.pep.xml
```

Or if both experiment file and identification file are in the same folder:

```
clustersheep --name=mysession /path/to/the/folder
```

If you have many experiment and identification files, you can compose a list of paths and save it as a text file, then load the file list with `--file-list=` option, like below:

```
clustersheep --name=mysession --file-list=/path/to/the/list.txt
```

where `list.txt` has following content:

```
/path/to/file-1.mzML
/path/to/file-2.pep.xml
/path/to/file-3.mzML
/path/to/file-4.pep.xml
          ⋮
```

<a name="usage-cluster-viewer"></a>
### Cluster viewer

After the finish of a clustering session, you can explore clusters in cluster viewer. To enter cluster viewer, execute ClusterSheep with `--stay-interactive` option, like below:

```
clustersheep --load-session=/path/to/mysession.cssess --stay-interactive
```

for loading an existing session or

```
clustersheep --name=mysession --stay-interactive path-to-files
```

for creating a new session.

<a name="usage-cluster-viewer-cluster-visualization"></a>
#### Visualization

Under cluster viewer, you can visualize a cluster by inputting its id. The cluster is drawn using [force-directed graph drawing](https://en.wikipedia.org/wiki/Force-directed_graph_drawing) algorithm. Under this algorithm, spectra that form more edges to its neighbours stay closed to each other while spectra that form less edges repel each other. The algorithm makes sub-clusters can be visually identified, for example, the cluster below.

![mixed-cluster](https://github.com/kpto/ClusterSheep/raw/release/docs/mixed-cluster.png)

In the above drawing, dots are nodes/spectra. An edge between two nodes means they have a similarity score higher than the threshold. If identifications are imported, identified spectra are colored. Spectra with the same identification are colored with the same colour. Black means no identification.

<a name="usage-cluster-viewer-spectrum-plotting"></a>
#### Spectrum plotting

When a cluster is drawn, you can interact with spectra by hovering your mouse pointer over a node, then press the following keys on your keyboard for different command:

* `T` to tag a spectrum.
* `D` to draw tagged spectrum/spectrum pair.
* `C` to clear tags.
* `V` to turn on or off verificative plotting.
* `I` to turn on or off identifications swapping of tagged spectrum pair.

For example, to plot a spectrum in a cluster, hover your mouse point over a node, press `T` button on your keyboard to tag the spectrum, then press `D` button to draw. The plot of the picked spectrum is shown on a new window. You can tag two spectra to plot them against each other like the image below. By default, spectra are plotted verificatively, meaning that the plotted spectrum is post-processed with mz range clipping, precursor peak removal and rank-transformation with the same configuration value used in clustering. If you want to plot the raw spectrum recorded in MS experiment files, press `V` to turn off verificative plotting or press again to turn it on again.

If the plotted spectrum is identified, the peptide sequence and the probability will be printed in the graph. Also, peaks of fragment ions are highlighted, as shown below. The theoretical spectrum referenced is generated using pyopenms. If two spectra are tagged and identifications swapping is turned on by pressing `I`, their identifications will be swapped and peaks are highlighted using the swapped identification.

![spectrum-plot](https://github.com/kpto/ClusterSheep/raw/release/docs/spectrum-plot.png)

Spectra cannot be plotted if input MS experiment files are moved because spectra data is read from source files.

<a name="usage-cluster-viewer-export"></a>
#### Export

Clusters can be exported by executing command `export` in cluster viewer. Information of nodes and edges will be written in a `Tab` delimited text file. If only a specific cluster is wanted, input the cluster ID (1000 for example) as argument like below:

```
export 1000
```

Or if you want to export all clusters, input argument `all`:

```
export all
```

File name can be specified by using `file=` argument:

```
export all file=path-to-destination
```

Before a cluster is exported, a process called `enrichment` is performed on that cluster first. This process append identification of all spectra in that cluster by searching the identification database to include such information in exported clusters. When exporting all clusters, such process will take a considerable time.

The exported file is sectioned. Section headers have a prefix of `#` followed by the section name. Each section is explained below.

```
# Columns // Titles of columns, matching the Tab delimited columns in Clusters section
    # Cluster
    # Nodes
    # Edges

# Clusters // Begin of exported clusters
    # Cluster
        Cluster ID and meta data, as shown in columns section
    # Nodes
        List of spectra
    # Edges
        List of edges, ID is the internal ID which corresponds to the first column
        of nodes section, NOT scan number
        
    # Cluster
    # Nodes
    # Edges
    
        ⋮
```

The following is an actual example:

```
# Columns
# Cluster
ID	Num of nodes	Num of edges	Num of identifications	Major identification	Identified ratio	Average precursor mass
# Nodes
ID	File	Scan num	Identification	Probability
# Edges
ID of source	ID of target	Dot product

# Clusters
# Cluster
5	4	6	1	n[33]SGK[160]VDVINAAK[160]/3	1.0	399.9366149902344
# Nodes
2950	/path/to/ms-experiment-file-1.mzXML	2378	n[33]SGK[160]VDVINAAK[160]/3	0.998657
2951	/path/to/ms-experiment-file-1.mzXML	2446	n[33]SGK[160]VDVINAAK[160]/3	0.99567
2953	/path/to/ms-experiment-file-2.mzXML	2758	n[33]SGK[160]VDVINAAK[160]/3	0.999604
2952	/path/to/ms-experiment-file-2.mzXML	2825	n[33]SGK[160]VDVINAAK[160]/3	0.999494
# Edges
2950	2951	0.71452552
2950	2952	0.72903901
2950	2953	0.75207931
2951	2952	0.77206755
2951	2953	0.77758884
2952	2953	0.82543987

# Cluster
7	3	2	0	None	0.0	402.2897033691406
# Nodes
3076	/path/to/ms-experiment-file-2.mzXML	7818	None	None
3079	/path/to/ms-experiment-file-3.mzXML	7886	None	None
3077	/path/to/ms-experiment-file-3.mzXML	7932	None	None
# Edges
3076	3077	0.81118226
3077	3079	0.71993017

    ⋮
```

Since export is done in parallel by multi-processing, clusters may not be written in order of cluster ID.

<a name="usage-advanced"></a>
### Advanced

<a name="usage-advanced-configuration"></a>
#### Configuration

There is a set of default values of clustering related parameters. For example, the precursor mass difference tolerance which is `1.1` and dot product threshold which is `0.7`. A comparison involving two spectra with their precursor mass difference larger than the tolerance is skipped as they are assumed to have different identity. An edge of two spectra is recorded only if the similarity score is higher than the threshold. To override these parameters, you can provide a configuration file like below:

```
clustersheep --name=mysession --config=/path/to/the/config-file.txt --file-list=/path/to/the/list.txt
```

where `config-file.txt` looks like this:

```
cg_precursor_tolerance = 2.0
cg_dot_product_threshold = 0.6
```

The above content changes the precursor mass difference to `2.0` and the dot product threshold to `0.6`. For more available parameters, generates a configuration template file by executing ClusterSheep with no parameter, as shown in section [Usage](#usage). The configuration template includes a detail explanation of all parameters.

Since ClusterSheep is built with traceable, the configuration used is recorded in a clustering session. You can use option `--print-session` to print the parameter values of it, like below:

```
clustersheep --load-session=/path/to/mysession.cssess --print-session
```

<a name="usage-advanced-material-preparation"></a>
#### Material preparation

In the case that you need to perform the same experiment multiple times but with different configuration, you may want to avoid the duplicated generations of index file `csindx` and rank-transformed spectra file `csrksp` to save time and disk space usage. To do so, you can run a new clustering session just to prepare the intermediate materials by using the option `--preparation-only`, like below:

```
clustersheep --name=prepsession --preparation-only --file-list=/path/to/the/list.txt
```

The above command will create an unfinished session. If you load this session, ClusterSheep will continue the flow and perform clustering on the prepared materials and finish the session. But the unfinished session is more useful if you fork it instead using option `--fork=`, like below:

```
clustersheep --load-session=/path/to/prepsession.cssess --fork=forkedsession --config=/path/to/config.txt
```

The above command will create a new session named by the value supplied by option `--fork=` (in the above example, `forkedsession`). The created session inherits everything from the loaded session including the prepared materials (index, rank-transformed spectra, logging and etc) by creating symbolic links to those files. By forking the preparation session with different configuration file, one can reuse the prepared material and run clustering on it with different configuration.

<a name="usage-advanced-python-interactive-console"></a>
#### Python interactive console

If you want to use Python to post-process a clustering session, you can turn on developer mode and enter Python interactive console from cluster viewer. To do so, load the session with option `--dev-mode` first:

```
clustersheep --load-session=/path/to/mysession.cssess --stay-interactive --dev-mode
```

You will enter cluster viewer as usual, but within developer mode, you can type `python` to enter Python interactive console. The prompt will be changed from `viewer >>>` to `python >>>`. All Python commands executed in the console will be logged into `cslogg` file as well.

<a name="usage-cluster-refinement"></a>
## Cluster refinement

Cluster refinement is a demonstrative feature to show the value of cluster structure. By computing [betweenness centrality](https://en.wikipedia.org/wiki/Betweenness_centrality) of each node, it identifies bridges within a cluster and remove them to produce clearer clusters. Using the same example cluster [shown above](#usage-cluster-viewer-cluster-visualization), two sub-clusters having different peptide identifications are bridged by the middle node. You can expect the betweenness centrality of that node will be a lot higher than other nodes and thus the node is detected and removed, releasing the two sub-clusters.

Since this feature is not officially a part of ClusterSheep, by default it is turned off. You can turn it on in a configuration file with parameter `cr_outlier_threshold`, see section [Configuration](#usage-advanced-configuration) for detail about configuration file.

<a name="limitation"></a>
## Limitation

Although ClusterSheep is built with big data in mind, the graph building process after pairwise similarity computation may still take unhandleable amount of ram if a cluster contains too many spectra. This is because ClusterSheep depends on graph-tool for all graph related processing and graph-tool does not support memory-mapping. Depending on a third party package reduces the development time of ClusterSheep significantly, but ClusterSheep will need a customized graph processing component that guarantees low memory use in the future. Before this happens, ClusterSheep is not truly a software that can handle big data.

If the value of dot product threshold is too low, all spectra will be clustered together and form a gigantic cluster. As long as the value is reasonable, clustering millions of spectra should not be a problem.

<a name="future"></a>
## Future

* Remove dependency of graph-tool.
* Overlap the clustering process and preparation process to hide the overhead of material preparation.
* GPU accelerated rank-transformation.
* Write a GUI.
* Test big data support.
* Decouple components.
* Write unit test.
* Supports OpenCL and AMD GPUs.
* Progressive clustering, add spectra after a clustering session.

<a name="bugs-or-requests"></a>
## Bugs or requests

For any bug reporting and feature request, send an email to [kpto@connect.ust.hk](mailto:kpto@connect.ust.hk).

<a name="thesis"></a>
## Thesis

[Building peptide spectral library by spectral clustering using graphics processing units (GPUs)](https://lbezone.ust.hk/bib/991012551465103412)