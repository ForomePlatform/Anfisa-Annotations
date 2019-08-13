package org.forome.annotation.connector.gnomad;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.forome.annotation.connector.gnomad.old.GnomadDataConnectorOld;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;

public class GnomadDataConnectorDiffTest {

    private final static Logger log = LoggerFactory.getLogger(GnomadDataConnectorDiffTest.class);

    @Test
    public void test() throws Exception {
        Assert.assertEquals("-CATCATCATCAT", GnomadDataConnectorOld.diff("CCATCATCATCAT", "C"));
        Assert.assertEquals("-CAT", GnomadDataConnectorOld.diff("CCATCAT", "CCAT"));
        Assert.assertEquals("-CATCATCAT", GnomadDataConnectorOld.diff("CCATCATCATCAT", "CCAT"));
        Assert.assertEquals("-CAT", GnomadDataConnectorOld.diff("CCATCAT", "CCAT"));
        Assert.assertEquals("-CATCAT", GnomadDataConnectorOld.diff("CCATCATCATCAT", "CCATCAT"));
        Assert.assertEquals("CAT", GnomadDataConnectorOld.diff("CCATCATCATCAT", "CCATCATCATCATCAT"));
    }

    @Test
    public void test1() throws Exception {
        Path file = Paths.get("/home/kris/processtech/tmp/bch0004_wgs_2.vcf");
        VCFFileReader vcfFileReader = new VCFFileReader(file, false);

        Iterator<VariantContext> iterator = vcfFileReader.iterator();
        while (iterator.hasNext()) {
            VariantContext variantContext = iterator.next();
            log.debug("Variant: {}", variantContext);
        }
        vcfFileReader.close();
    }

}
