## Guidelines for working with Ensembl VEP locally using its docker image.

1. Copy either env-GRCh37.sh or env-GRCh38.sh into env.sh

2. Edit `env.sh` file and define all parameters for the run. The parameters that may change between different builds: 

		VEP_RELEASE=release_105.0       # Ensemble VEP release
		VEP_ASSEMBLY=GRCh38             # The homo-sapiens species genome assembly
        VEP_PORT=5306                   # vep docker image specific port

3. Launch build script - this should be done only once for a particular release and genome, when it is required to download the new cache and fasta files. They will be saved into a local directory, defined by `$VEP_CACHE_DIR` variable (`/data/vep/vep-cache-${VEP_ASSEMBLY}-${VEP_RELEASE}` by default).

		bash build.sh

4. Now it is possible to run VEP, providing the input file as an argument:

                bash run.sh pgp3140_wgs_hlpanel.vcf

    An example that shows how to run it in detached mode from an arbitrary location:

                nohup bash ~/work/Anfisa-Annotations/pipeline/projects/ensembl-vep/run.sh pgp3140_wgs_hlpanel.vcf &> vep.log &

    In case it is needed to switch to a different release/assembly - edit `env.sh` again and launch `run.sh`.
