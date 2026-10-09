# Linux/NVIDIA runtime adapted from the pinned RFdiffusion upstream Dockerfile.
# Build: docker build --platform linux/amd64 -f containers/Dockerfile \
#          -t binderflow-rfdiffusion:86507b6-cu116 .
FROM nvcr.io/nvidia/cuda:11.6.2-cudnn8-runtime-ubuntu20.04

ARG RFDIFFUSION_REV=86507b6538f51fce57b5a72477165f03999ed7ae
LABEL org.opencontainers.image.source="https://github.com/RosettaCommons/RFdiffusion"
LABEL org.opencontainers.image.revision="${RFDIFFUSION_REV}"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      git ca-certificates python3.9 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && git init /app/RFdiffusion \
    && cd /app/RFdiffusion \
    && git remote add origin https://github.com/RosettaCommons/RFdiffusion.git \
    && git fetch --depth 1 origin ${RFDIFFUSION_REV} \
    && git checkout --detach FETCH_HEAD

RUN python3.9 -m pip install --no-cache-dir pip==24.3.1 \
    && python3.9 -m pip install --no-cache-dir \
      numpy==1.23.5 scipy==1.10.1 \
      dgl==1.0.2+cu116 -f https://data.dgl.ai/wheels/cu116/repo.html \
      torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116 \
      e3nn==0.3.3 wandb==0.12.0 protobuf==3.20.3 pynvml==11.0.0 decorator==5.1.0 \
      hydra-core==1.3.2 pyrsistent==0.19.3 \
      git+https://github.com/NVIDIA/dllogger@0478734ff7be75adde8d160e04872664d1c62e5f \
      /app/RFdiffusion/env/SE3Transformer \
    && python3.9 -m pip install --no-cache-dir --no-deps /app/RFdiffusion \
    && python3.9 -c "import torch, dgl, rfdiffusion"

ENV DGLBACKEND=pytorch
WORKDIR /app/RFdiffusion
# Nextflow supplies the command; weights are staged from --models at runtime.
