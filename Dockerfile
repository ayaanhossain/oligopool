# Use an Ubuntu LTS base so Python package availability stays stable.
FROM ubuntu:24.04

# Avoid interactive prompts during package installs.
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and minimal build tooling (some dependencies may require compilation).
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set up a virtual environment for all Python installs.
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv "$VIRTUAL_ENV"
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Set docker build directory.
WORKDIR /oligopool-docker

# Copy the library files into the container.
COPY . /oligopool-docker

# Install the library (and ensure console scripts like `oligopool` / `op` are available).
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir .

# Install Jupyter in the same environment.
RUN pip install --no-cache-dir jupyter

# Install additional utilities useful in interactive containers.
RUN apt-get update && apt-get install -y --no-install-recommends \
    nano \
    vim \
    tree \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Change prompt anchor (interactive shells source ~/.bashrc).
RUN echo 'export PS1="\[\033[0;34m\][\[\033[1;32m\]\u\[\033[0;34m\]]-(\[\033[1;34m\]\w\[\033[0;34m\])\n \[\033[1;36m\]>> \[\033[0m\]"' >> /root/.bashrc

# Default working directory (mounted by docs/docker-notes.md examples).
WORKDIR /op-workspace

# Expose ports for Jupyter (default: 8080).
EXPOSE 8080 8081 8082

# Start an interactive bash shell by default.
ENTRYPOINT ["/bin/bash"]
