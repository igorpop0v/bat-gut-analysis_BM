FROM ubuntu:24.04

RUN apt-get update && \
    apt-get install -y wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/sortmerna && \
    wget \
    https://github.com/sortmerna/sortmerna/releases/download/v7.0.0/sortmerna-7.0.0-Linux.sh \
    -O /tmp/sortmerna-install.sh && \
    bash /tmp/sortmerna-install.sh \
        --skip-license \
        --prefix=/opt/sortmerna \
        --exclude-subdir && \
    rm /tmp/sortmerna-install.sh

ENV PATH="/opt/sortmerna/bin:${PATH}"

ENTRYPOINT ["sortmerna"]
