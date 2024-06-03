
from rocker/shiny-verse:4.3.3
run apt-get update && \
  # for igraph
  apt-get install -y libxml2-dev libglpk-dev

run R -e "install.packages(c('dplyr', 'reshape2', 'RColorBrewer', 'pROC', 'text2vec', 'flexdashboard'))"

run apt-get update && apt-get install -y git
run R -e "remotes::install_git('https://gitlab.com/thomaschln/sgraph.git', ref = 'c512896fa646b6be56568834600f473272e4e53c')"


add ./ /kgraph
run R -e "devtools::install('kgraph', dependencies = TRUE)"
