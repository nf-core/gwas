# MPH 0.55.1: upstream static Linux amd64 binary and corresponding source.
# Debian supplies /bin/bash; procps supports Nextflow task tracing.
ARG BASE_IMAGE=docker.io/library/debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171

FROM ${BASE_IMAGE} AS fetch

ARG TARGETARCH
ARG MPH_VERSION=0.55.1
ARG MPH_COMMIT=13ffe633e0170d35e01e6164fc2bb4e01f35e2f2

ARG MPH_ZIP_SHA256=9a2d976baa0e281b097ed17d20d622111000a8d0fa040157a5b7540c16bd4454
ARG MPH_BIN_SHA256=02b083e6c2ed9100e13236386faa45138f5d6a0ed570c994b504af72484a60eb
ARG MPH_SRC_SHA256=4c49ea8d8e08a7e7b132c129650bbafae54e79b8aedd3acf4106b89f1b4a7ad3
ARG MKL_LICENSE_WHL_SHA256=1f9514cec32ab1e467c7fc44674485622d9b151e91b9318e4d97249a28595671
ARG IOMP_LICENSE_WHL_SHA256=adfb32b0dde6b3a95ced62608b05a45a7cfad928054136ddbfe0563e16cbf33c

RUN set -eux; \
    arch="${TARGETARCH:-amd64}"; \
    if [ "${arch}" != "amd64" ]; then \
        echo "ERROR: MPH ${MPH_VERSION} is distributed by upstream only as a" >&2; \
        echo "       statically linked Intel MKL binary for Linux x86_64." >&2; \
        echo "       Refusing to build an untested image for TARGETARCH=${arch}." >&2; \
        echo "       Build with --platform=linux/amd64." >&2; \
        exit 1; \
    fi

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates curl unzip; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN set -eux; \
    curl -fsSL -o mph.zip \
      "https://github.com/jiang18/mph/releases/download/v${MPH_VERSION}/mph-v${MPH_VERSION}-x86_64-linux.zip"; \
    echo "${MPH_ZIP_SHA256}  mph.zip" | sha256sum -c -; \
    unzip -q -j mph.zip mph -d .; \
    echo "${MPH_BIN_SHA256}  mph" | sha256sum -c -; \
    chmod 0755 mph; \
    head -c 20 mph | od -An -tx1 | tr -d ' \n' \
      | grep -q '^7f454c4602010103000000000000000002003e00'; \
    rm -f mph.zip

RUN set -eux; \
    curl -fsSL -o mph-source.tar.gz \
      "https://github.com/jiang18/mph/archive/refs/tags/v${MPH_VERSION}.tar.gz"; \
    echo "${MPH_SRC_SHA256}  mph-source.tar.gz" | sha256sum -c -; \
    mkdir -p /licenses/source; \
    tar -xzf mph-source.tar.gz -C /licenses/source --strip-components=1; \
    cp /licenses/source/LICENSE /licenses/GPL-3.0-only.txt; \
    grep -q 'GNU GENERAL PUBLIC LICENSE' /licenses/GPL-3.0-only.txt; \
    grep -q 'Version 3, 29 June 2007' /licenses/GPL-3.0-only.txt; \
    mv mph-source.tar.gz /licenses/mph-v${MPH_VERSION}-source.tar.gz

RUN set -eux; \
    curl -fsSL -o mkl-include.whl \
      "https://files.pythonhosted.org/packages/eb/ac/076148c4edd7dadda49ee7504bd4a8319e125a675bc28f3966244e3e3812/mkl_include-2023.0.0-py2.py3-none-manylinux1_i686.whl"; \
    echo "${MKL_LICENSE_WHL_SHA256}  mkl-include.whl" | sha256sum -c -; \
    curl -fsSL -o intel-openmp.whl \
      "https://files.pythonhosted.org/packages/26/49/ee0d58c2e73e6faee43c2e45b271e68e1079bdc49a44f847e83d2afae134/intel_openmp-2023.0.0-py2.py3-none-manylinux1_i686.whl"; \
    echo "${IOMP_LICENSE_WHL_SHA256}  intel-openmp.whl" | sha256sum -c -; \
    unzip -q -o -j mkl-include.whl \
      'mkl_include-2023.0.0.dist-info/LICENSE.txt' -d /tmp/mkl; \
    unzip -q -o -j intel-openmp.whl \
      'intel_openmp-2023.0.0.dist-info/LICENSE.txt' \
      'intel_openmp-2023.0.0.data/data/licensing/compiler/openmp/third-party-programs.txt' -d /tmp/iomp; \
    mv /tmp/mkl/LICENSE.txt  /licenses/Intel-Simplified-Software-License-Oct2022.oneMKL.txt; \
    mv /tmp/iomp/LICENSE.txt /licenses/Intel-EULA-Developer-Tools-Oct2021.libiomp5.txt; \
    mv /tmp/iomp/third-party-programs.txt /licenses/Intel-OpenMP-third-party-programs.txt; \
    grep -q 'Intel Simplified Software License (Version October 2022)' \
      /licenses/Intel-Simplified-Software-License-Oct2022.oneMKL.txt; \
    grep -q 'Intel End User License Agreement for Developer Tools' \
      /licenses/Intel-EULA-Developer-Tools-Oct2021.libiomp5.txt; \
    rm -f mkl-include.whl intel-openmp.whl

RUN set -eux; \
    find /licenses -type d -exec chmod 0755 {} +; \
    find /licenses -type f -exec chmod 0644 {} +

FROM ${BASE_IMAGE}

ARG MPH_VERSION=0.55.1
ARG MPH_COMMIT=13ffe633e0170d35e01e6164fc2bb4e01f35e2f2

LABEL org.opencontainers.image.title="mph" \
      org.opencontainers.image.description="MPH ${MPH_VERSION} - MINQUE/REML partitioning of heritability from multiple genomic relationship matrices (upstream prebuilt static Intel MKL x86_64 binary)" \
      org.opencontainers.image.version="${MPH_VERSION}" \
      org.opencontainers.image.revision="${MPH_COMMIT}" \
      org.opencontainers.image.source="https://github.com/jiang18/mph" \
      org.opencontainers.image.url="https://jiang18.github.io/mph/" \
      org.opencontainers.image.documentation="https://jiang18.github.io/mph/" \
      org.opencontainers.image.licenses="GPL-3.0-only AND LicenseRef-IntelSimplifiedSoftwareOct2022 AND LicenseRef-IntelEULADeveloperToolsOct2021" \
      org.opencontainers.image.vendor="Jicai Jiang (upstream); repackaged for nf-core/gwas" \
      org.opencontainers.image.base.name="docker.io/library/debian:bookworm-slim" \
      org.opencontainers.image.base.digest="sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171" \
      io.nfcore.arch="linux/amd64 only (static Intel MKL build; requires SSE4.2)"

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends procps; \
    rm -rf /var/lib/apt/lists/*

COPY --from=fetch --chmod=0755 /build/mph /usr/local/bin/mph
COPY --from=fetch /licenses/ /usr/share/doc/mph/

RUN set -eux; \
    mph > /tmp/mph-banner.txt 2>&1 || true; \
    grep -q 'MPH by Jicai Jiang' /tmp/mph-banner.txt; \
    grep -q "Version ${MPH_VERSION} " /tmp/mph-banner.txt; \
    test -s /usr/share/doc/mph/GPL-3.0-only.txt; \
    test -s /usr/share/doc/mph/Intel-Simplified-Software-License-Oct2022.oneMKL.txt; \
    test -d /usr/share/doc/mph/source/src; \
    rm -f /tmp/mph-banner.txt

WORKDIR /data
CMD ["mph"]
