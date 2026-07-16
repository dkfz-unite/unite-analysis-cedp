FROM ghcr.io/dkfz-unite/docker-rdev-biocmanager:latest AS base

FROM base AS install
COPY ./src/install /src
WORKDIR /src
RUN apt-get update && \
apt-get install -y --no-install-recommends libuv1-dev && \
rm -rf /var/lib/apt/lists/* && \
Rscript install.R
RUN apt-get clean

FROM install AS final
COPY ./src/run /src
COPY ./app /app
WORKDIR /app
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
ENV ASPNETCORE_hostBuilder:reloadConfigOnChange=false
ENV UNITE_COMMAND="Rscript"
ENV UNITE_COMMAND_ARGUMENTS="run.R {data}/{proc}"
ENV UNITE_SOURCE_PATH="/src"
ENV UNITE_DATA_PATH="/mnt/data"
ENV UNITE_PROCESS_LIMIT="1"
EXPOSE 80
CMD ["/app/commands", "--urls", "http://0.0.0.0:80"]
