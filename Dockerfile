FROM ubuntu:20.04 as builder

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# ensure any pipes fail correctly, may impact other things
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# work in an area that will get discarded
WORKDIR /tmp/build

# Prevent interactive options during installs
ENV DEBIAN_FRONTEND=noninteractive

# apt stuff
# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends software-properties-common \
&& add-apt-repository ppa:deadsnakes/ppa \
&& apt-get -yq install --no-install-recommends python3.9 \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3.9 1

# not in final image
# hadolint ignore=DL3008
RUN apt-get -yq install --no-install-recommends \
    build-essential \
    python3.9-dev \
    python3.9-distutils \
    python3.9-venv \
    curl \
&& curl -sSL --retry 10 https://bootstrap.pypa.io/get-pip.py > get-pip.py \
&& python3.9 get-pip.py

# base environment
ENV OPT /opt/wsi-t78
ENV VIRTUAL_ENV=$OPT/venv
RUN python3.9 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY README.md .
COPY pygas/ pygas/
COPY setup.py requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY tests/ tests/
COPY examples/ examples/
RUN python3.9 ./setup.py develop && tests/scripts/run_unit_tests.sh

# we the coverage report (and junit.xml) in for info and later use
RUN mkdir -p /var/www/pygas && mv htmlcov /var/www/pygas/test-coverage && mv junit.xml /var/www/pygas/.

# as we use COPY we can cleanup stuff from the testing layers and actually get a smaller image
RUN pip uninstall -y pre-commit pytest-cov pytest coverage

# deploy properly
RUN python3.9 ./setup.py install

########################## FINAL IMAGE ##########################
FROM ubuntu:20.04

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Prevent interactive options during installs
ENV DEBIAN_FRONTEND=noninteractive

# apt stuff
# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends software-properties-common \
&& add-apt-repository ppa:deadsnakes/ppa \
&& apt-get -yq install --no-install-recommends python3.9 \
# make sure all security patches are applied
&& apt-get -yq update && apt-get install -yq --no-install-recommends unattended-upgrades \
&& unattended-upgrade -d -v \
&& apt-get remove -yq unattended-upgrades \
&& apt-get autoremove -yq \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/* \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3.9 1

ENV OPT /opt/wsi-t78
ENV VIRTUAL_ENV=$OPT/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY --from=builder $OPT $OPT
# the test coverate html output - we can grab and pass to gitlab-pages
COPY --from=builder /var/www/ /var/www/

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
