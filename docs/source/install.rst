Installation
============

To install this workflow, you have to clone the github repository as described below. This workflow uses Snakemake and Conda, you can install them following their manual instructions as linked in the sections below or using the docker image provided by the developpers team of Snakemake.

After, the installation of Snakemake and Conda, you have to install our wrapper system following the command described below. The wrapper system allows to facilitate the re-usability of programs and to simplify the use in Snakemake. All wrappers are available in the `API section <https://translatome.github.io/TRAIN/end_toc.html>`_.


Clone workflow into desired working directory
---------------------------------------------

Copy the Translatome workflow into the desired folder *path/to/workdir*:

.. code-block:: shell

    git clone https://github.com/Translatome/Translatome.git path/to/workdir


Installation of Snakemake through MiniConda
-------------------------------------------

For the installation of Conda, please refer to the `official documentation <https://conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_ and for the installation of Mamba, please also refer to the `official documentation <https://mamba.readthedocs.io/en/latest/installation.html>`_.

Here, this is an example:

1. Download the installer miniconda: `https://conda.io/miniconda.html <https://conda.io/miniconda.html>`_

2. In your Terminal window, run:

.. code-block:: shell

    bash Miniconda3-latest-Linux-x86_64.sh

3. Follow the prompts on the installer screens.

4. Install Snakemake in a dedicated environment

.. code-block:: shell

    conda install -c bioconda -c conda-forge -n snakemake snakemake python=3.10
    conda activate snakemake


.. Note::
    Conda is an open source package management system and environment management system that runs on Windows, macOS, and Linux. It allows the users to easily install and manage installation of programs. The wanted programs are listed in YAML environment files. A conda environment is a directory containing the collection of wanted programs, it doesn't require to be the administrator of your computer or cluster for the installation. Snakemake can be installed through Conda.


Installation of Snakemake through Docker
----------------------------------------

The official Docker image can be found at: `https://hub.docker.com/r/snakemake/snakemake <https://hub.docker.com/r/snakemake/snakemake>`_

This image contains Conda, Mamba and Snakemake.

For more information on Docker images, read the `official documentation <https://docs.docker.com/get-started/overview/>`_.

.. Note::
    Docker is an open platform for developing, shipping, and running applications; it provides the ability to package and run an application in a loosely isolated environment called a container. Docker containers can run on a developer's local laptop, on physical or virtual machines in a data center, on cloud providers, or in a mixture of environments. Containers with pre-installed programs are available on the Docker Hub as "images". The user can find an image containing Snakemake and Conda.


Installation of wrapper system package
--------------------------------------

In order to install the wrappers on your computer, go to the *translatome/* folder and use the:

.. code-block:: shell
    
    # require pip for python 3
    pip install .

.. Note::
    This command line install the wrapper locally, if any modification is made to the wrappers, you have to re-install them using the same command line.

