#!/bin/bash

vcf=$1
header=${vcf}.header
body=${vcf}.body
fixed="${vcf%.*}"_fixed.vcf

grep -E "^#" ${vcf} > ${header}
grep -vE "^#" ${vcf} > ${body}
cp ${header} ${fixed}
awk '$5!="."' ${body} >> ${fixed}
bgzip ${fixed}
tabix ${fixed}.gz

rm {body}
rm {header}


