# BUILD:
# git clone --depth 1 -b sltp-lite --single-branch git@github.com:aig-upf/fs-private.git fs
# git clone --depth 1 -b d2l --single-branch git@github.com:aig-upf/sltp.git sltp
# git clone --depth 1 -b devel --single-branch git@github.com:aig-upf/tarski.git tarski
# cd fs && git submodule update --init && cd..
# docker build -t sltp .
#
# docker save sltp | bzip2 | pv | ssh awscluster 'bunzip2 | docker load'
#
# RUN:
#    docker run -it sltp  # Run interactive shell
#    docker run -it sltp /root/projects/sltp/experiments/d2l/run.py blocks:clear

FROM ubuntu:focal
MAINTAINER Guillem Francès guillem.frances@upf.edu


ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/London"
RUN apt-get update && apt-get install --no-install-recommends -y \
	build-essential \
	cmake \
	ca-certificates \
	curl \
	python3 python3-pip python3-dev python3-setuptools \
	scons \
	git \
	pkg-config \
	libboost-program-options-dev libboost-filesystem-dev libboost-system-dev \
	libboost-chrono-dev libboost-timer-dev libboost-serialization-dev \
    # OpenWbo:
    vim libz-dev libgmp-dev \
    && rm -rf /var/lib/apt/lists/*


# Set up some environment variables
ENV HOME /root
ENV BASE_DIR "${HOME}/projects"
ENV CXX g++

# Install the wheel module for subsequents pip installs
RUN pip3 install wheel

RUN mkdir -p $BASE_DIR

#################################
# Install & build OpenWbo
#################################
ENV OPENWBO_GITHUB_URL "https://github.com/sat-group/open-wbo/archive/master.tar.gz"
RUN cd $BASE_DIR && curl -SL $OPENWBO_GITHUB_URL | tar -xz \
	&& mv open-wbo-master open-wbo \
	&& cd open-wbo \
	&& make rs \
    && cp ./open-wbo_static /bin/ \
    && make clean

#################################
# Install & build the FS planner
#################################
ENV FS_PATH "${BASE_DIR}/fs"
ADD fs $FS_PATH
RUN cd $FS_PATH && git submodule update --init && python3 ./build.py -p


#################################
# Install & build Tarski
#################################
ENV TARSKI_PATH "${BASE_DIR}/tarski"
ADD tarski $TARSKI_PATH
RUN cd $TARSKI_PATH && pip3 install -e .


#################################
# Install & build SLTP
#################################
ENV SLTP_PATH "${BASE_DIR}/sltp"
ADD sltp $SLTP_PATH
RUN cd $SLTP_PATH && pip3 install -e .
RUN cd $SLTP_PATH/src/generators && cmake . && make -j4 && rm -rf CMakeFiles


WORKDIR $SLTP_PATH
#ENTRYPOINT ["sh", "-c", "$SLTP_PATH/experiments/run.py"]
#CMD ["bash"]

COPY ./docker-entrypoint.sh /
ENTRYPOINT ["/docker-entrypoint.sh"]
