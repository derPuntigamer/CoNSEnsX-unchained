FROM ubuntu:18.10 AS builder
LABEL maintainer="dudola.daniel@itk.ppke.hu"

ENV LANG=C.UTF-8

RUN apt-get update -q \
    && apt-get install -y -q wget python3 python3-pip wget \
    && wget users.itk.ppke.hu/~dudda/progs_consensx.tar.gz \
    && pip3 install pipenv

ENV PYROOT /pyroot
ENV PYTHONUSERBASE $PYROOT

WORKDIR /build

COPY Pipfile Pipfile.lock ./

RUN PIP_USER=1 PYTHONUSERBASE=$PYROOT pipenv install --system --deploy

#                           .ÄÄÄ-.
#                          (___/\ \
#        ,                 (|^ ^ ) )    At this point, the dependencies 
#       /(                _)_\=_/  (    are built in the builder container.
# ,..__/ `\          ____(_/_ ` \   )   
#  `\    _/        _/---._/(_)_  `\ (   Copy dependencies to final image.
#    '--\ `-.__..-'    /.    (_), |  )
#        `._        ___\_____.'_| |__/
#           `~----"`   `-.........'
FROM ubuntu:18.10

ENV PGPASSWORD=password
ENV LANG=C.UTF-8

COPY --from=builder /progs_consensx.tar.gz .

# add i386 archutecture (pales dependency) and install software used by CoNSEnsX
RUN dpkg --add-architecture i386 \
    && apt-get update -q \
    && apt-get install -y -q postgresql-client libc6:i386 libncurses5:i386 \
                             libstdc++6:i386 libx11-6:i386 python3 python3-setuptools \
    && mkdir -p /usr/src/app \
    && tar -C /usr/src/app -zxf progs_consensx.tar.gz \
    && rm -rf progs_consensx.tar.gz \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

ENV PYROOT /pyroot
ENV PYTHONPATH $PYROOT/lib/python:$PATH
ENV PYTHONUSERBASE $PYROOT

COPY --from=builder $PYROOT/lib/ $PYROOT/lib/ 
COPY . /usr/src/app/consensx.itk.ppke.hu

WORKDIR /usr/src/app/consensx.itk.ppke.hu

CMD ["bash", "docker_cmd.sh"]

EXPOSE 8000
