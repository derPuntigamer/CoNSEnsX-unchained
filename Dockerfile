FROM ubuntu:18.10
LABEL maintainer="dudola.daniel@itk.ppke.hu"

ENV DEBIAN_FRONTEND=noninteractive
ENV PGPASSWORD=password
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# add i386 archutecture (pales dependency) and install software used by CoNSEnsX
RUN dpkg --add-architecture i386 \
    && apt-get update -q \
    && apt-get install -y -q postgresql-client libc6:i386 libncurses5:i386 libstdc++6:i386 libx11-6:i386 wget python3 python3-distutils python3-setuptools python3-pip \
    && pip3 install pipenv \
    && mkdir -p /usr/src/app \
    && wget users.itk.ppke.hu/~dudda/progs_consensx.tar.gz \
    && tar -C /usr/src/app -zxvf progs_consensx.tar.gz \
    && rm progs_consensx.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

COPY . /usr/src/app/consensx.itk.ppke.hu
WORKDIR /usr/src/app/consensx.itk.ppke.hu

RUN pipenv install

CMD ["bash", "docker_cmd.sh"]

EXPOSE 8000