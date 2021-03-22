# Start with neurodebian trusty 16.04

FROM neurodebian:xenial
MAINTAINER Hua Wu <huawu@stanford.edu>

# Install dependencies
RUN echo deb http://neurodeb.pirsquared.org data main contrib non-free >> /etc/apt/sources.list.d/neurodebian.sources.list
RUN echo deb http://neurodeb.pirsquared.org xenial main contrib non-free >> /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-get update && apt-get -y install \
    python python-pip \
    afni \
    fsl-core

RUN pip install --upgrade pip==20.3 \
    && pip install numpy==1.16.4 \
    && pip install nibabel==2.1.0 
#    && pip install nipype==0.10.0

# Make directory for flywheel spec
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}
COPY run.py ${FLYWHEEL}/run
COPY manifest.json ${FLYWHEEL}/manifest.json

# Put script into flywheel folder
COPY tsnr.py ${FLYWHEEL}/tsnr.py

# Set the entrypoint
ENTRYPOINT ["/flywheel/v0/run"]
