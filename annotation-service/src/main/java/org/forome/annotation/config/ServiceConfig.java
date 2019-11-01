package org.forome.annotation.config;

import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import org.forome.annotation.config.connector.*;
import org.forome.annotation.config.ensemblvep.EnsemblVepConfig;
import org.forome.annotation.config.frontend.FrontendConfig;
import org.forome.annotation.config.notification.NotificationSlackConfig;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public class ServiceConfig {

    public final FrontendConfig frontendConfig;

    public final EnsemblVepConfig ensemblVepConfigConnector;

    public final GnomadConfigConnector gnomadConfigConnector;
    public final ClinVarConfigConnector clinVarConfigConnector;
    public final HgmdConfigConnector hgmdConfigConnector;
    public final GTFConfigConnector gtfConfigConnector;
    public final SpliceAIConfigConnector spliceAIConfigConnector;
    public final ConservationConfigConnector conservationConfigConnector;
    public final RefConfigConnector refConfigConnector;
    public final GTEXConfigConnector gtexConfigConnector;

    public final NotificationSlackConfig notificationSlackConfig;

    public ServiceConfig() throws Exception {
        this(Paths.get("config.json").toAbsolutePath());
    }

    public ServiceConfig(Path configFile) throws Exception {
        if (!Files.exists(configFile)) {
            throw new RuntimeException("File: " + configFile.toString() + " not found");
        }
        JSONObject configFileJson;
        try (InputStream is = Files.newInputStream(configFile, StandardOpenOption.READ)) {
            configFileJson = (JSONObject) new JSONParser(JSONParser.DEFAULT_PERMISSIVE_MODE).parse(is);
        }

        frontendConfig = new FrontendConfig((JSONObject) configFileJson.get("frontend"));

        ensemblVepConfigConnector = new EnsemblVepConfig((JSONObject) configFileJson.get("ensembl-vep"));

        JSONObject jConnectors = (JSONObject) configFileJson.get("connectors");
        gnomadConfigConnector = new GnomadConfigConnector((JSONObject) jConnectors.get("gnomad"));
        clinVarConfigConnector = new ClinVarConfigConnector((JSONObject) jConnectors.get("clinvar"));
        hgmdConfigConnector = new HgmdConfigConnector((JSONObject) jConnectors.get("hgmd"));
        gtfConfigConnector = new GTFConfigConnector((JSONObject) jConnectors.get("gtf"));
        spliceAIConfigConnector = new SpliceAIConfigConnector((JSONObject) jConnectors.get("spliceai"));
        conservationConfigConnector = new ConservationConfigConnector((JSONObject) jConnectors.get("conservation"));
        refConfigConnector = new RefConfigConnector((JSONObject) jConnectors.get("ref"));
        gtexConfigConnector = new GTEXConfigConnector((JSONObject) jConnectors.get("gtex"));

        JSONObject jNotifications = (JSONObject) configFileJson.get("notification");
        if (jNotifications != null) {
            JSONObject jNotificationSlackConfig = (JSONObject) jNotifications.get("slack");
            if (jNotificationSlackConfig != null) {
                notificationSlackConfig = new NotificationSlackConfig(jNotificationSlackConfig);
            } else {
                notificationSlackConfig = null;
            }
        } else {
            notificationSlackConfig = null;
        }
    }

}
