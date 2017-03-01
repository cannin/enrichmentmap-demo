FROM cannin/r-shiny-server:ubuntu-14.04.4_r-3.3.2_java-8_shiny-server-1.5.2.837
MAINTAINER cannin

COPY r-requirements.txt /
COPY installPackages.R /
COPY runInstallPackages.R /
RUN R -e 'source("runInstallPackages.R")'

# Copy sample apps
COPY * /srv/shiny-server/emdemo/

# Expose Shiny server
EXPOSE 3838
#CMD ["shiny-server"]
