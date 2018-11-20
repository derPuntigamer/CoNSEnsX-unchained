FROM ubuntu:18.04
LABEL maintainer="dudola.daniel@itk.ppke.hu"

ENV DEBIAN_FRONTEND=noninteractive
ENV PGPASSWORD=password

# add i386 archutecture (pales dependency) and install software used by CoNSEnsX
RUN dpkg --add-architecture i386 \
    && apt-get update && apt-get install -y \
    git postgresql-client libc6:i386 libncurses5:i386 libstdc++6:i386 libx11-6:i386 \
    wget python3 python3-django python3-numpy \
    python3-scipy python3-matplotlib python3-dev python3-psycopg2 \
    && rm -rf /var/lib/apt/lists/*

# create directory for CoNSEnsX and cd to it
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

# download and install ProDy
RUN wget https://pypi.python.org/packages/1e/9d/1afe15ccbadc847305ea92044c0011d77349e89a61cf34b4b3b3d60ff28d/ProDy-1.9.tar.gz#md5=be40dd72d89a97f0fd0fe1739b93578f \
    && tar -zxvf ProDy-1.9.tar.gz
WORKDIR /usr/src/app/ProDy-1.9
RUN python3 setup.py build && python3 setup.py install

WORKDIR /usr/src/app

RUN wget users.itk.ppke.hu/~dudda/progs_consensx.tar.gz \
    && tar -zxvf progs_consensx.tar.gz

# add this and below command will run without cache
# this is because you don't have to reinstall everything upon a GIT push event
ARG CACHEBUST=1

# let's get down to business
# RUN git clone -b master https://github.com/PPKE-Bioinf/consensx.itk.ppke.hu
COPY . /usr/src/app/consensx.itk.ppke.hu
WORKDIR /usr/src/app/consensx.itk.ppke.hu

CMD ["bash", "docker_cmd.sh"]

EXPOSE 8000