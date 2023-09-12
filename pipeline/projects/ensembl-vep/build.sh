#!/bin/bash

## Run the initial procedure to download and cache Ensembl VEP data for the particular release

set -euxo pipefail

VEP_SETUP_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd )

source ${VEP_SETUP_DIR}/env.sh

if [ -d $VEP_CACHE_DIR ]; then
    echo "Forome VEP setup is already built. Remove $VEP_CACHE_DIR cache directory to generate a new build."
    exit 1
fi      

echo "Pull Ensembl VEP image for ${VEP_RELEASE} release"

docker pull $VEP_IMAGE

echo "Download VEP cache files"

mkdir -p $VEP_CACHE_DIR
chmod a+rwx $VEP_CACHE_DIR

docker run --rm -it --init --user $(id -u):$(id -g) -e 'HOME=/opt/vep' -v $(realpath $VEP_CACHE_DIR):/opt/vep/.vep $VEP_IMAGE perl INSTALL.pl -a p -g $VEP_PLUGINS
docker run --rm -it --init --user $(id -u):$(id -g) -e 'HOME=/opt/vep' -v $(realpath $VEP_CACHE_DIR):/opt/vep/.vep $VEP_IMAGE perl INSTALL.pl -a cf -s $VEP_SPECIES -y $VEP_ASSEMBLY

echo "Download Plugins config files"

VEP_RELEASE_SHORT=$(echo ${VEP_RELEASE}|awk -F'[_.]' '{print$2}')

wget https://raw.github.com/Ensembl/VEP_plugins/release/${VEP_RELEASE_SHORT}/ExACpLI_values.txt -o /dev/null -O ${VEP_CACHE_DIR}/Plugins/ExACpLI_values.txt
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz -o /dev/null && tar xzf fordownload.tar.gz -C ${VEP_CACHE_DIR}/Plugins/ && rm fordownload.tar.gz
wget https://raw.github.com/Ensembl/VEP_plugins/release/${VEP_RELEASE_SHORT}/LoFtool_scores.txt -o /dev/null -O ${VEP_CACHE_DIR}/Plugins/LoFtool_scores.txt

echo "Forome VEP setup is done"
