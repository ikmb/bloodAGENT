LABEL description="Bloodagent development image"

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
    libbz2-1.0 \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Arbeitsverzeichnis setzen
WORKDIR /app

# Build-Output vom vorherigen Schritt kopieren
COPY --from=builder /app/dist/Release/GNU-Linux/bloodAGENT /app/bloodAGENT
COPY --from=builder /app/external/htslib/libhts.so /app/libhts.so
RUN ln -s /app/libhts.so /app/libhts.so.3
RUN mkdir /licenses
COPY Third_Party_Licenses.md /licenses/
COPY LICENSE /licenses/
COPY --from=builder /app/external/libBigWig/libBigWig.so /app/libBigWig.so

# Setze LD_LIBRARY_PATH, damit die Shared Libraries gefunden werden
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/


# Standardkommando setzen
CMD ["./bloodAGENT"]
