package org.forome.annotation.struct.variant.cnv;

import org.forome.annotation.struct.Chromosome;
import org.forome.annotation.struct.variant.VariantType;
import org.forome.annotation.struct.variant.vep.VariantVep;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class VariantCNV extends VariantVep {

    public final Map<String, Genotype> genotypes;

    public VariantCNV(Chromosome chromosome, int start, int end, List<Genotype> genotypes) {
        super(chromosome, start, end);
        this.genotypes = Collections.unmodifiableMap(
                genotypes.stream()
                        .collect(Collectors.toMap(item -> item.sampleName, item -> item))
        );
    }

    public Genotype getGenotype(String sample) {
        return genotypes.get(sample);
    }

    @Override
    public VariantType getVariantType() {
        return VariantType.CNV;
    }

    @Override
    public String getRef() {
        return "DEL";
    }

    @Override
    public List<String> getAltAllele() {
        return Collections.singletonList("-");
    }
}
