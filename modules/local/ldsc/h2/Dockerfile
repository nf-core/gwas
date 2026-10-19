# syntax=docker/dockerfile:1.7

FROM python:3.9.13-slim-bullseye@sha256:dcf2eafca55558d8b1aa73edd6aa41b7187c5bcb63e533a7b04a0673f81f37fe AS build

ARG LDSC_COMMIT="6c673952cee74bd5c57aef1555a03b1c015399a0"
ARG LDSC_ARCHIVE_SHA256="9b3f0c26c181fdcd5c668183fa5ef6f4d0df28a371831f49f1b6254402ce2321"

LABEL org.opencontainers.image.authors="nf-core/gwas contributors" \
      org.opencontainers.image.description="CBIIT Python 3 LD Score Regression command-line runtime" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.revision="${LDSC_COMMIT}" \
      org.opencontainers.image.source="https://github.com/CBIIT/ldsc" \
      org.opencontainers.image.version="3.0.2"

ADD --checksum=sha256:${LDSC_ARCHIVE_SHA256} \
    https://github.com/CBIIT/ldsc/archive/${LDSC_COMMIT}.tar.gz \
    /tmp/ldsc.tar.gz

RUN apt-get update \
    && apt-get install --yes --no-install-recommends build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /tmp/ldsc \
    && tar --extract --gzip --file /tmp/ldsc.tar.gz --directory /tmp/ldsc --strip-components=1 \
    && python -m venv /opt/venv \
    && /opt/venv/bin/python -m pip install --no-cache-dir \
        bitarray==2.6.0 \
        nose==1.3.7 \
        numpy==1.23.3 \
        pandas==1.5.0 \
        pybedtools==0.9.1 \
        pysam==0.19.1 \
        python-dateutil==2.8.2 \
        pytz==2022.4 \
        scipy==1.9.2 \
        six==1.16.0 \
    && /opt/venv/bin/python -m pip install --no-cache-dir --no-deps /tmp/ldsc \
    && sed -n "s/^__version__ = ['\"]\\([^'\"]*\\)['\"]/\\1/p" /tmp/ldsc/ldsc.py > /opt/venv/ldsc-native-version \
    && printf '%s\n' "${LDSC_COMMIT}" > /opt/venv/ldsc-source-revision

FROM python:3.9.13-slim-bullseye@sha256:dcf2eafca55558d8b1aa73edd6aa41b7187c5bcb63e533a7b04a0673f81f37fe

ARG LDSC_COMMIT="6c673952cee74bd5c57aef1555a03b1c015399a0"

LABEL org.opencontainers.image.authors="nf-core/gwas contributors" \
      org.opencontainers.image.description="CBIIT Python 3 LD Score Regression command-line runtime" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.revision="${LDSC_COMMIT}" \
      org.opencontainers.image.source="https://github.com/CBIIT/ldsc" \
      org.opencontainers.image.version="3.0.2"

COPY --from=build /opt/venv /opt/venv

ENV PATH="/opt/venv/bin:${PATH}"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/* \
    && ldsc.py --help >/dev/null \
    && munge_sumstats.py --help >/dev/null \
    && test "$(python -c 'import importlib.metadata; print(importlib.metadata.version("ldsc"))')" = "3.0.2" \
    && test "$(cat /opt/venv/ldsc-native-version)" = "3.0.1" \
    && test "$(cat /opt/venv/ldsc-source-revision)" = "${LDSC_COMMIT}"

CMD ["/bin/bash"]
