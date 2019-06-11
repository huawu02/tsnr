# Start with neurodebian trusty 14.04

FROM neurodebian:trusty
MAINTAINER Hua Wu <huawu@stanford.edu>

# Install dependencies
RUN echo deb http://neurodeb.pirsquared.org data main contrib non-free >> /etc/apt/sources.list.d/neurodebian.sources.list
RUN echo deb http://neurodeb.pirsquared.org trusty main contrib non-free >> /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-get update && apt-get -y install \
    python-dev python-pip \
    fsl-5.0-core
RUN pip install --upgrade pip \
    && pip install numpy==1.11.0 \
    && pip install nibabel==2.1.0
#Run pip install nipype==0.10.0

# Make directory for flywheel spec
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}
COPY run.py ${FLYWHEEL}/run
COPY manifest.json ${FLYWHEEL}/manifest.json

# Put script into flywheel folder
COPY tsnr.py ${FLYWHEEL}/tsnr.py

# Set the entrypoint
ENTRYPOINT ["/flywheel/v0/run"]
