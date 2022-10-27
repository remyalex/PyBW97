FROM ubuntu:16.04

MAINTAINER Remy Galan "remyalexander@gmail.com"

RUN apt-get update -y && \
    apt-get install -y qgis python-qgis

COPY . /pybw

WORKDIR /pybw

#ENTRYPOINT [ "python" ]

#CMD [ "./code/mmiMesh.py" ]
