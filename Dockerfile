FROM jupyter/datascience-notebook:python-3.8.8


USER root

COPY ./vafator /home/jovyan/vafator/
COPY ./setup.py /home/jovyan
COPY ./setup.cfg /home/jovyan
COPY ./requirements.txt /home/jovyan
COPY ./MANIFEST.in /home/jovyan
COPY ./README.md /home/jovyan
COPY ./notebooks /home/jovyan/notebooks/

# install system dependencies
RUN apt-get update && apt-get install -y zlib1g-dev

# install vafator
RUN python -m pip install --upgrade pip
RUN python setup.py bdist_wheel
RUN pip install dist/*

# installs some extensions to Jupyter
RUN jupyter labextension install jupyterlab-plotly --no-build \
   && jupyter labextension install @jupyterlab/toc-extension \
   && jupyter labextension install @jupyterlab/git \
   && jupyter lab build -y \
   && jupyter lab clean -y

