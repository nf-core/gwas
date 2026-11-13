# METASOFT 2.0.1 with its native jar and bundled RE2 p-value table.
# Debian provides /bin/bash; procps supports Nextflow task tracing.
FROM docker.io/library/debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171

ARG METASOFT_VERSION=2.0.1
ARG METASOFT_URL=https://hanlab-snu.github.io/METASOFT/meta/repository/2.0.1/Metasoft.zip
ARG METASOFT_SHA256=56ccaf5e04c3a7c4f51a98c44f495ee512399f73057504146928d0b1d7d148df
ARG METASOFT_SRC_URL=https://hanlab-snu.github.io/METASOFT/repository/2.0.0/Metasoft_src.zip
ARG METASOFT_SRC_SHA256=60b05b1f3d4a38eff3488456af546d81b5a928169010bed35b47d7b2002d2625

LABEL org.opencontainers.image.title="METASOFT" \
      org.opencontainers.image.version="2.0.1" \
      org.opencontainers.image.description="METASOFT 2.0.1 — fixed-effects, DerSimonian-Laird random-effects, Han & Eskin RE2 and binary-effects meta-analysis for GWAS" \
      org.opencontainers.image.authors="Buhm Han" \
      org.opencontainers.image.url="https://hanlab-snu.github.io/METASOFT/index.html" \
      org.opencontainers.image.source="https://hanlab-snu.github.io/METASOFT/repository/2.0.0/Metasoft_src.zip" \
      org.opencontainers.image.licenses="LicenseRef-METASOFT-Buhm-Han AND LicenseRef-CERN-COLT" \
      org.metasoft.license.file="/usr/share/doc/metasoft/LICENSE" \
      org.metasoft.redistribution.permission="granted-by-author-email"

COPY metasoft /usr/local/bin/metasoft

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends openjdk-17-jre-headless procps ca-certificates curl unzip; \
    rm -rf /var/lib/apt/lists/*; \
    mkdir -p /opt/metasoft /usr/share/doc/metasoft; \
    cd /tmp; \
    curl -fsSL -o Metasoft.zip "${METASOFT_URL}"; \
    curl -fsSL -o Metasoft_src.zip "${METASOFT_SRC_URL}"; \
    echo "${METASOFT_SHA256}  Metasoft.zip"         | sha256sum -c -; \
    echo "${METASOFT_SRC_SHA256}  Metasoft_src.zip" | sha256sum -c -; \
    unzip -q Metasoft.zip -d unpacked; \
    install -m 0444 unpacked/Metasoft.jar              /opt/metasoft/Metasoft.jar; \
    install -m 0444 unpacked/HanEskinPvalueTable.txt   /opt/metasoft/HanEskinPvalueTable.txt; \
    install -m 0444 unpacked/example.txt               /opt/metasoft/example.txt; \
    install -m 0444 unpacked/LICENSE                   /usr/share/doc/metasoft/LICENSE; \
    install -m 0444 unpacked/README                    /usr/share/doc/metasoft/README; \
    # CERN COLT's own notice, shipped inside the fat jar's One-JAR wrapper docs \
    unzip -q -o -j unpacked/Metasoft.jar 'doc/one-jar-license.txt' -d /usr/share/doc/metasoft; \
    # Retain the pristine source archive so the copyright notice in the .java \
    # headers travels with the binary, as the licence requires. \
    install -m 0444 Metasoft_src.zip /usr/share/doc/metasoft/Metasoft_src-2.0.0.zip; \
    printf '%s\n' "${METASOFT_VERSION}" > /opt/metasoft/VERSION; \
    chmod 0444 /opt/metasoft/VERSION; \
    chmod 0555 /usr/local/bin/metasoft; \
    rm -rf /tmp/Metasoft.zip /tmp/Metasoft_src.zip /tmp/unpacked; \
    # Fail the build unless the baked artefacts run offline and reproduce \
    # a known example result. \
    cd /tmp; \
    metasoft -input /opt/metasoft/example.txt -output /tmp/selftest.out -log /tmp/selftest.log; \
    grep -q '^rs9	3	1.03788E-48' /tmp/selftest.out; \
    rm -f /tmp/selftest.out /tmp/selftest.log

# NB: deliberately no JAVA_TOOL_OPTIONS — setting it, even to "", makes the JVM
# print "Picked up JAVA_TOOL_OPTIONS:" on stderr on every invocation.
ENV METASOFT_HOME=/opt/metasoft \
    METASOFT_JAR=/opt/metasoft/Metasoft.jar \
    METASOFT_PVALUE_TABLE=/opt/metasoft/HanEskinPvalueTable.txt

WORKDIR /work
CMD ["metasoft", "-help"]
