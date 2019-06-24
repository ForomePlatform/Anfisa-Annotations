package org.forome.annotation.connector.conservation.struct;

import com.google.common.annotations.VisibleForTesting;

import java.util.Locale;

public class Conservation {

    public final Double priPhCons;
    public final Double mamPhCons;
    public final Double verPhCons;
    public final Double priPhyloP;
    public final Double mamPhyloP;
    public final Double verPhyloP;
    public final Double gerpRS;
    private final Double gerpRSpval;
    public final Double gerpN;
    public final Double gerpS;

    public Conservation(
            Double priPhCons, Double mamPhCons,
            Double verPhCons, Double priPhyloP,
            Double mamPhyloP, Double verPhyloP,
            Double gerpRS, Double gerpRSpval,
            Double gerpN, Double gerpS
    ) {
        this.priPhCons = priPhCons;
        this.mamPhCons = mamPhCons;
        this.verPhCons = verPhCons;
        this.priPhyloP = priPhyloP;
        this.mamPhyloP = mamPhyloP;
        this.verPhyloP = verPhyloP;
        this.gerpRS = gerpRS;
        this.gerpRSpval = gerpRSpval;
        this.gerpN = gerpN;
        this.gerpS = gerpS;
    }

    public String getGerpRSpval(){
        if (gerpRSpval == null) {
            return null;
        }
        return convFromL(gerpRSpval);
    }


    @VisibleForTesting
    public static Double convToL(String value){
        if (value == null) {
            return null;
        }
        if ("0".equals(value)) {
            return -10000D;
        }
        String[] splits = value.split("e", 2);
        if (!splits[1].startsWith("-")) {
            throw new NumberFormatException();
        }
        return Math.log10(Double.parseDouble(splits[0])) + Integer.parseInt(splits[1]);
    }

    @VisibleForTesting
    public static String convFromL(Double value){
        if (value == null) {
            return null;
        }
        if (value < -9999) {
            return "0";
        }
        double power = Math.floor(value);
        return String.format(Locale.ENGLISH, "%.5fe%s", Math.pow(10D, value-power), (long)power);
    }
}

