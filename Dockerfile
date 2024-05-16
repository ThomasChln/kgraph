
from rocker/shiny-verse:4.3.3
run apt-get update && \
  # for igraph
  apt-get install -y libxml2-dev libglpk-dev

run R -e "install.packages(c('dplyr', 'reshape2', 'RColorBrewer', 'pROC', 'text2vec', 'flexdashboard'))"

run apt-get install -y git
run R -e "remotes::install_git('https://gitlab.com/thomaschln/sigmagraph.git', ref = 'f05000a4e2e11f38f984f70e824abb461f281699')"

add ./ /kgraph
run R -e "devtools::install('kgraph', dependencies = TRUE)"
