# Basis-Image mit Build-Tools
FROM ubuntu:22.04 AS builder

# Abhängigkeiten installieren
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    git \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Arbeitsverzeichnis setzen
WORKDIR /app

# Quellcode kopieren
COPY . .

# Kompilieren (ersetze 'main.cpp' mit deinem eigentlichen Code)
RUN cd ./external/htslib && make clean  && make && \
    cd ../libBigWig && make clean  && make && \
    cd ../..  && make clean && make CONF=Release

# Runtime-Image (kleineres Image)
FROM ubuntu:22.04

# Laufzeit-Abhängigkeiten installieren
RUN apt-get update && apt-get install -y \
    zlib1g \
    liblzma5 \
    libhts3 \
    libbz2-1.0 \
    libcurl4 \
    libbigwig0 \
    && rm -rf /var/lib/apt/lists/*

# Arbeitsverzeichnis setzen
WORKDIR /app

# Build-Output vom vorherigen Schritt kopieren
COPY --from=builder /app/dist/Release/GNU-Linux/deepblood /app/bloodAGENT

# Setze LD_LIBRARY_PATH, damit die Shared Libraries gefunden werden
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/external/htslib:/app/external/libBigWig


# Standardkommando setzen
CMD ["./bloodAGENT"]
