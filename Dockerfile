# EMLaunch - Electromagnetic Launcher Simulation
# Dockerfile for running Jupyter notebooks with full environment

FROM ubuntu:24.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    python3 \
    python3-pip \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Julia 1.12 from official tarball
ARG JULIA_VERSION=1.12.5
RUN ARCH=$(dpkg --print-architecture) && \
    if [ "$ARCH" = "arm64" ]; then JULIA_ARCH="aarch64"; JULIA_URL_ARCH="aarch64"; \
    else JULIA_ARCH="x86_64"; JULIA_URL_ARCH="x64"; fi && \
    curl -fsSL "https://julialang-s3.julialang.org/bin/linux/${JULIA_URL_ARCH}/1.12/julia-${JULIA_VERSION}-linux-${JULIA_ARCH}.tar.gz" \
    | tar -xz -C /opt && \
    ln -s /opt/julia-${JULIA_VERSION}/bin/julia /usr/local/bin/julia

# Install Jupyter
RUN pip3 install --break-system-packages --no-cache-dir \
    jupyter \
    notebook

# Install IJulia kernel
RUN julia -e 'using Pkg; Pkg.add("IJulia")'

# Set working directory
WORKDIR /app

# Copy project files
COPY Project.toml ./
COPY src/ ./src/
COPY notebooks/ ./notebooks/
COPY LICENSE README.md ./

# Install Julia dependencies
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Expose Jupyter port
EXPOSE 8888

# Set working directory to notebooks
WORKDIR /app/notebooks

# Start Jupyter notebook server
CMD ["jupyter", "notebook", \
     "--ip=0.0.0.0", \
     "--port=8888", \
     "--no-browser", \
     "--allow-root", \
     "--NotebookApp.token=''", \
     "--NotebookApp.password=''"]
