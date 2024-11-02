# Use Ubuntu as the base image
FROM ubuntu:latest

# Avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3.12 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Set up a virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Set docker directory
WORKDIR /oligopool-docker

# Copy the library files into the container
COPY . /oligopool-docker

# Install the library and its dependencies
RUN pip install --no-cache-dir .

# Install any additional packages you want available in the environment
RUN apt-get update && apt-get install -y \
    nano \
    vim \
    tree \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]