package org.forome.annotation.struct.variant;

import org.forome.annotation.struct.Chromosome;

import java.util.List;

public abstract class Variant {

    public final Chromosome chromosome;
    public final int start;
    public final int end;

    public Variant(Chromosome chromosome, int start, int end) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    public abstract VariantType getVariantType();

    public abstract Genotype getGenotype(String sample);

    public abstract String getRef();

    public abstract List<String> getAltAllele();

    public abstract String getMostSevereConsequence();

}
