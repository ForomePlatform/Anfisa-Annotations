package org.forome.annotation.service.cnv;

import net.minidev.json.JSONObject;
import org.forome.annotation.Main;
import org.forome.annotation.annotator.input.CNVFileIterator;
import org.forome.annotation.config.ServiceConfig;
import org.forome.annotation.connector.anfisa.AnfisaConnector;
import org.forome.annotation.connector.anfisa.struct.AnfisaInput;
import org.forome.annotation.connector.anfisa.struct.AnfisaResult;
import org.forome.annotation.connector.clinvar.ClinvarConnector;
import org.forome.annotation.connector.conservation.ConservationConnector;
import org.forome.annotation.connector.gnomad.GnomadConnector;
import org.forome.annotation.connector.gnomad.GnomadConnectorImpl;
import org.forome.annotation.connector.gtf.GTFConnector;
import org.forome.annotation.connector.hgmd.HgmdConnector;
import org.forome.annotation.connector.liftover.LiftoverConnector;
import org.forome.annotation.connector.spliceai.SpliceAIConnector;
import org.forome.annotation.service.ensemblvep.EnsemblVepService;
import org.forome.annotation.service.ensemblvep.external.EnsemblVepExternalService;
import org.forome.annotation.service.ssh.SSHConnectService;
import org.forome.annotation.struct.variant.cnv.VariantCNV;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.nio.file.Paths;

public class CNVMain {

    private final static Logger log = LoggerFactory.getLogger(Main.class);

    public static void main(String[] args) throws Exception {
        ServiceConfig serviceConfig = new ServiceConfig();
        SSHConnectService sshTunnelService = new SSHConnectService();
        GnomadConnector gnomadConnector = new GnomadConnectorImpl(sshTunnelService, serviceConfig.gnomadConfigConnector, (t, e) -> crash(e));
        SpliceAIConnector spliceAIConnector = new SpliceAIConnector(sshTunnelService, serviceConfig.spliceAIConfigConnector, (t, e) -> crash(e));
        ConservationConnector conservationConnector = new ConservationConnector(sshTunnelService, serviceConfig.conservationConfigConnector);
        HgmdConnector hgmdConnector = new HgmdConnector(sshTunnelService, serviceConfig.hgmdConfigConnector);
        ClinvarConnector clinvarConnector = new ClinvarConnector(sshTunnelService, serviceConfig.clinVarConfigConnector);
        LiftoverConnector liftoverConnector = new LiftoverConnector();
        GTFConnector gtfConnector = new GTFConnector(sshTunnelService, serviceConfig.gtfConfigConnector, (t, e) -> crash(e));
        EnsemblVepService ensemblVepService = new EnsemblVepExternalService((t, e) -> crash(e));
        AnfisaConnector anfisaConnector = new AnfisaConnector(
                gnomadConnector,
                spliceAIConnector,
                conservationConnector,
                hgmdConnector,
                clinvarConnector,
                liftoverConnector,
                gtfConnector,
                (t, e) -> crash(e)
        );

        Path pathVcf = Paths.get("/home/kris/processtech/tmp/3/cnv.vcf");
        CNVFileIterator cnvFileIterator = new CNVFileIterator(pathVcf);

        while (cnvFileIterator.hasNext()) {
            VariantCNV variant = cnvFileIterator.next();
            JSONObject vepJson = ensemblVepService.getVepJson(variant, "-").get();
            AnfisaInput anfisaInput = new AnfisaInput.Builder().build();
            AnfisaResult anfisaResult = anfisaConnector.build(
                    null, anfisaInput, variant, vepJson );
            log.debug("anfisaResult: " + anfisaResult);
        }

        log.debug("end");
    }

    public static void crash(Throwable e) {
        log.error("Application crashing ", e);
        System.exit(1);
    }
}
