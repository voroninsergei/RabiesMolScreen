FROM mambaorg/micromamba:1.5.3

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n rabiesmol -f /tmp/environment.yml && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "rabiesmol", "/bin/bash", "-c"]

WORKDIR /app
COPY . /app

CMD ["/bin/bash"]
