/*
 Copyright (c) 2019. Vladimir Ulitin, Partners Healthcare and members of Forome Association

 Developed by Vladimir Ulitin and Michael Bouzinier

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

	 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

package org.forome.annotation.config;

import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import org.forome.annotation.config.connector.*;
import org.forome.annotation.config.database.DatabaseConfig;
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

	public final DatabaseConfig databaseConfig;

	public final EnsemblVepConfig ensemblVepConfigConnector;

	public final GnomadConfigConnector gnomadConfigConnector;
	public final ClinVarConfigConnector clinVarConfigConnector;
	public final HgmdConfigConnector hgmdConfigConnector;
	public final GTFConfigConnector gtfConfigConnector;
	public final SpliceAIConfigConnector spliceAIConfigConnector;
	public final RefConfigConnector refConfigConnector;
	public final GTEXConfigConnector gtexConfigConnector;
	public final PharmGKBConfigConnector pharmGKBConfigConnector;

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

		databaseConfig = new DatabaseConfig((JSONObject) configFileJson.get("database"));

		ensemblVepConfigConnector = new EnsemblVepConfig((JSONObject) configFileJson.get("ensembl-vep"));

		JSONObject jConnectors = (JSONObject) configFileJson.get("connectors");
		gnomadConfigConnector = new GnomadConfigConnector((JSONObject) jConnectors.get("gnomad"));
		clinVarConfigConnector = new ClinVarConfigConnector((JSONObject) jConnectors.get("clinvar"));
		hgmdConfigConnector = new HgmdConfigConnector((JSONObject) jConnectors.get("hgmd"));
		gtfConfigConnector = new GTFConfigConnector((JSONObject) jConnectors.get("gtf"));
		spliceAIConfigConnector = new SpliceAIConfigConnector((JSONObject) jConnectors.get("spliceai"));
		refConfigConnector = new RefConfigConnector((JSONObject) jConnectors.get("ref"));
		gtexConfigConnector = new GTEXConfigConnector((JSONObject) jConnectors.get("gtex"));
		pharmGKBConfigConnector = new PharmGKBConfigConnector((JSONObject) jConnectors.get("pharmgkb"));

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
