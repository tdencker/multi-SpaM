FROM alpine:latest as multispam-builder
RUN apk add --no-cache gcc g++ cmake ninja

WORKDIR /opt/multispam
ADD . .

WORKDIR /opt/multispam/build
RUN cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
RUN ninja

FROM alpine:latest
LABEL maintainer="Thomas Dencker <thomas.dencker@stud.uni-goettingen.de>"
RUN apk add --no-cache libstdc++ python libgomp

COPY --from=multispam-builder /opt/multispam/bin/multi-SpaM /usr/local/bin/multi-SpaM
COPY --from=multispam-builder /opt/multispam/bin/max-cut-tree /usr/local/bin/max-cut-tree
COPY --from=multispam-builder /opt/multispam/multispam.py /usr/local/multispam.py

# TODO: wrapper for avx/sse3

ENTRYPOINT ["python", "/usr/local/multispam.py"]