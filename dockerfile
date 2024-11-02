# Use Ubuntu as the base image
FROM ubuntu:latest

# Avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.12 -- update later as necessary
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
    htop \
    && rm -rf /var/lib/apt/lists/*

# Make nice CLI prompt
RUN echo 'export PS1="\[\033[0;34m\][\[\033[1;32m\]\u\[\033[0;34m\]]─(\[\033[1;34m\]\w\[\033[0;34m\])\n \[\033[1;36m\]>> \[\033[0m\]"' >> ~/.bashrc

# Set working directory
WORKDIR /op-workspace

# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]