# Debian supplies /bin/bash directly; procps supports Nextflow task tracing.
FROM docker.io/library/debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171

ARG LDAK_VERSION="6.3"
ARG LDAK_REVISION="b755ab71b6dcf8f36f40bb64c799ec858d2ca492"
ARG LDAK_SHA256="8f8f5102acf1c001641322b17cd3dd624c1710c645e8a673235f6a250f11d0c6"

LABEL org.opencontainers.image.authors="Andong Ni <lyh970817@users.noreply.github.com>"
LABEL org.opencontainers.image.description="Nextflow-compatible image containing the official LDAK 6.3 Linux binary"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"
LABEL org.opencontainers.image.source="https://github.com/dougspeed/LDAK"
LABEL org.opencontainers.image.revision="${LDAK_REVISION}"
LABEL org.opencontainers.image.version="${LDAK_VERSION}"
LABEL org.opencontainers.image.ldak.binary-sha256="${LDAK_SHA256}"

ADD --chmod=0755 --checksum=sha256:8f8f5102acf1c001641322b17cd3dd624c1710c645e8a673235f6a250f11d0c6 \
    https://raw.githubusercontent.com/dougspeed/LDAK/b755ab71b6dcf8f36f40bb64c799ec858d2ca492/ldak6.3.linux \
    /usr/local/bin/ldak

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN echo "${LDAK_SHA256}  /usr/local/bin/ldak" | sha256sum -c - \
    && mkdir -p /opt/ldak \
    && printf '%s\n' "${LDAK_VERSION}" > /opt/ldak/version \
    && printf '%s\n' "${LDAK_REVISION}" > /opt/ldak/revision \
    && printf '%s\n' "${LDAK_SHA256}" > /opt/ldak/sha256
