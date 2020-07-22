FROM lc_app_upper_base_04jul2020
WORKDIR /home/shiny
COPY webapp.R .
COPY tbl.twE12.deg.RData .
COPY twe12.bed .
COPY twtw.bed .
COPY pkgs/igvShiny_0.99.91.tar.gz .
RUN R CMD INSTALL igvShiny_0.99.91.tar.gz
COPY pkgs/ghdb_1.0.2.tar.gz .
RUN R CMD INSTALL ghdb_1.0.2.tar.gz
COPY pkgs/shinyModules_1.0.13.tar.gz .
RUN R CMD INSTALL shinyModules_1.0.13.tar.gz
COPY installPrerequisites.R .
RUN R -f installPrerequisites.R
# CMD ["R", "-e", "shiny::runApp('/home/shiny/lcApp.R', port=3838, host='0.0.0.0')"]
