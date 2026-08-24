# Dockerfile for MOCCa CI/CD on GitHub Actions
# 
# This image contains all build and runtime dependencies for MOCCa,
# eliminating the need for apt-get operations during workflow execution.
#
# Usage in workflows:
#   jobs:
#     compile:
#       runs-on: ubuntu-latest
#       container: ghcr.io/your-org/mocca-ci:latest

FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install all dependencies in one layer to minimize image size
# (Combining RUN commands reduces layers but increases layer size slightly)
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        gfortran make git \
        liblapack-dev libblas-dev \
        libhdf5-dev hdf5-tools \
        python3 python3-pip python3-numpy \
        time bc ca-certificates \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install specific scipy version via pip
RUN pip3 install --no-cache-dir numpy scipy==1.17.1

# Set up working directory
WORKDIR /workspace

# Verify installation
RUN gfortran --version && \
    python3 -c "import numpy; import scipy; print('scipy version:', scipy.__version__)" && \
    h5cc --version

# Label the image
LABEL maintainer="wouter.ryssens@ulb.be"
LABEL description="MOCCa CI/CD environment with all dependencies pre-installed"
