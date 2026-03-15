# bloodAGENT base image for Bioscope platform
# Builds bloodAGENT from source and bundles runtime dependencies.
# This image can be used standalone or extended by the genomics-blood-type-processor.

# Stage 1: Build bloodAGENT from source
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy source and initialize submodules
COPY . .
RUN git submodule update --init --recursive

RUN cd external/htslib && make && \
    cd ../libBigWig && make && \
    cd ../.. && make CONF=Release

# Stage 2: Runtime image with bloodAGENT binary + bioinformatics tools
FROM ubuntu:22.04
LABEL description="bloodAGENT base image for Bioscope genomics-blood-type-processor"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    liblzma5 \
    libbz2-1.0 \
    libcurl4 \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# bloodAGENT binary and shared libraries
COPY --from=builder /build/dist/Release/GNU-Linux/bloodAGENT /usr/local/bin/bloodAGENT
COPY --from=builder /build/external/htslib/libhts.so /usr/local/lib/libhts.so
COPY --from=builder /build/external/libBigWig/libBigWig.so /usr/local/lib/libBigWig.so
RUN ln -s /usr/local/lib/libhts.so /usr/local/lib/libhts.so.3 && ldconfig

# bloodAGENT configuration and annotation data
COPY --from=builder /build/data /data

# Licenses (BSD 2-Clause compliance)
COPY --from=builder /build/LICENSE /licenses/bloodAGENT-LICENSE
COPY --from=builder /build/Third_Party_Licenses.md /licenses/bloodAGENT-Third_Party_Licenses.md

# Default entrypoint for standalone usage; overridden by the processor image
ENTRYPOINT ["/usr/local/bin/bloodAGENT"]
CMD ["--help"]
