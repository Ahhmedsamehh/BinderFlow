FROM binderflow-rfdiffusion:86507b6-cu116
ARG PROTEINMPNN_REV=8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
RUN git init /app/ProteinMPNN \
    && cd /app/ProteinMPNN \
    && git remote add origin https://github.com/dauparas/ProteinMPNN.git \
    && git fetch --depth 1 origin ${PROTEINMPNN_REV} \
    && git checkout --detach FETCH_HEAD \
    && python3.9 -m pip install --no-cache-dir biopython==1.85 \
    && test -s vanilla_model_weights/v_48_020.pt
LABEL org.opencontainers.image.source="https://github.com/dauparas/ProteinMPNN"
LABEL org.opencontainers.image.revision="${PROTEINMPNN_REV}"
