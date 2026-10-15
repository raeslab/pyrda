# Dev image for pyrda: Python (uv-managed) + R (with RLdbRDA) + Claude Code.
# R base: rocker/r-ver (Ubuntu noble) ships with Posit Package Manager binaries,
# so CRAN packages like vegan install as prebuilt binaries.
FROM rocker/r-ver:4.5

COPY --from=ghcr.io/astral-sh/uv:0.12.10 /uv /uvx /usr/local/bin/

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

# System packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git curl ca-certificates sudo \
        libxml2 libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# R packages: RLdbRDA + its dependencies, plus languageserver for VS Code
# RLdbRDA is pinned to a commit for reproducible comparisons; bump deliberately.
ARG RLDBRDA_REF=98704d538d1ab5f186ebed5f192ab9febe3efd6c
RUN install2.r --error --skipinstalled --ncpus -1 \
        remotes vegan dplyr ggplot2 tibble tidyr magrittr languageserver \
    && Rscript -e "remotes::install_github('raeslab/RLdbRDA', ref = '${RLDBRDA_REF}', upgrade = 'never')" \
    && Rscript -e "library(RLdbRDA)" \
    && rm -rf /tmp/downloaded_packages

# Non-root user (Ubuntu noble ships a default 'ubuntu' user on UID 1000; remove it)
RUN (userdel -r ubuntu 2>/dev/null || true) \
    && groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME} \
    && echo "${USERNAME} ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME} \
    && mkdir -p /workspace && chown ${USERNAME}:${USERNAME} /workspace

USER ${USERNAME}

# Python via uv; venv lives outside the bind mount
ENV UV_LINK_MODE=copy \
    UV_PYTHON=3.13 \
    UV_PYTHON_PREFERENCE=only-managed \
    UV_PROJECT_ENVIRONMENT=/home/vscode/.venv \
    PATH="/home/vscode/.venv/bin:/home/vscode/.local/bin:$PATH"

RUN uv python install 3.13

# Claude Code (native installer)
RUN curl -fsSL https://claude.ai/install.sh | bash \
    && mkdir -p /home/vscode/.claude

WORKDIR /workspace
