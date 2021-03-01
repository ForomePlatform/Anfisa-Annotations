#!/bin/bash

## Run the initial procedure to download and cache Ensembl VEP data for the particular release

set -euxo pipefail

source env.sh

if [ -d $VEP_CACHE_DIR ]; then
    echo "Forome VEP setup is already built. Remove $VEP_CACHE_DIR cache directory to perform a new build."
    exit
fi      

echo "Pull Ensembl VEP image for ${VEP_RELEASE} release"

docker pull $VEP_IMAGE

echo "Download VEP cache files"

mkdir -p $VEP_CACHE_DIR
chmod a+rwx $VEP_CACHE_DIR

docker run --rm -it --init --user $(id -u):$(id -g) -e 'HOME=/opt/vep' -v $(realpath $VEP_CACHE_DIR):/opt/vep/.vep $VEP_IMAGE perl INSTALL.pl -a p -g $VEP_PLUGINS
docker run --rm -it --init --user $(id -u):$(id -g) -e 'HOME=/opt/vep' -v $(realpath $VEP_CACHE_DIR):/opt/vep/.vep $VEP_IMAGE perl INSTALL.pl -a cf -s $VEP_SPECIES -y $VEP_ASSEMBLY

echo "Download Plugins config files"

wget https://raw.github.com/Ensembl/VEP_plugins/${VEP_RELEASE/_/\/}/ExACpLI_values.txt -o /dev/null -O ${VEP_CACHE_DIR}/Plugins/ExACpLI_values.txt
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz -o /dev/null && tar xzf fordownload.tar.gz -C ${VEP_CACHE_DIR}/Plugins/ && rm fordownload.tar.gz
wget https://raw.github.com/Ensembl/VEP_plugins/${VEP_RELEASE/_/\/}/LoFtool_scores.txt -o /dev/null -O ${VEP_CACHE_DIR}/Plugins/LoFtool_scores.txt

#echo "export VEP_IMAGE=${VEP_IMAGE}   # VEP image" > ./_vep_env.sh
#echo "export VEP_CACHE_DIR=${VEP_CACHE_DIR}   # The local directory to store cached files" >>  ./_vep_env.sh

echo "Forome VEP setup is done"
