# EMLaunch - Electromagnetic Launcher Simulation
# Dockerfile for running Jupyter notebooks with full environment

FROM julia:1.12

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Jupyter
RUN pip3 install --break-system-packages --no-cache-dir \
    jupyter \
    notebook

# Install IJulia kernel
RUN julia -e 'using Pkg; Pkg.add("IJulia")'

# Set working directory
WORKDIR /app

# Copy project files
COPY Project.toml Manifest.toml ./
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
