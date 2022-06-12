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

package org.forome.annotation;

import org.forome.querypool.*;
import org.forome.annotation.config.Config;
import org.forome.annotation.config.ServiceConfig;
import org.forome.annotation.data.DatabaseConnector;
import org.forome.annotation.data.anfisa.AnfisaConnector;
import org.forome.annotation.data.clinvar.ClinvarConnector;
import org.forome.annotation.data.clinvar.mysql.ClinvarConnectorMysql;
import org.forome.annotation.data.gnomad.GnomadConnectorImpl;
import org.forome.annotation.data.gnomad.datasource.http.GnomadDataSourceHttp;
import org.forome.annotation.data.gtex.GTEXConnector;
import org.forome.annotation.data.gtex.mysql.GTEXConnectorMysql;
import org.forome.annotation.data.gtf.GTFConnector;
import org.forome.annotation.data.gtf.GTFConnectorImpl;
import org.forome.annotation.data.gtf.datasource.mysql.GTFDataConnector;
import org.forome.annotation.data.hgmd.HgmdConnector;
import org.forome.annotation.data.hgmd.mysql.HgmdConnectorMysql;
import org.forome.annotation.data.pharmgkb.PharmGKBConnector;
import org.forome.annotation.data.pharmgkb.mysql.PharmGKBConnectorMysql;
import org.forome.annotation.data.spliceai.SpliceAIConnector;
import org.forome.annotation.data.spliceai.SpliceAIConnectorImpl;
import org.forome.annotation.data.spliceai.datasource.http.SpliceAIDataSourceHttp;
import org.forome.annotation.database.DatabaseService;
import org.forome.annotation.database.entityobject.user.UserReadable;
import org.forome.annotation.exception.AnnotatorException;
import org.forome.annotation.exception.QueryPoolExceptionBuilder;
import org.forome.annotation.network.NetworkService;
import org.forome.annotation.network.component.UserEditableComponent;
import org.forome.annotation.service.database.DatabaseConnectService;
import org.forome.annotation.service.ensemblvep.EnsemblVepService;
import org.forome.annotation.service.ensemblvep.external.EnsemblVepExternalService;
import org.forome.annotation.service.notification.NotificationService;
import org.forome.annotation.service.source.SourceService;
import org.forome.annotation.service.ssh.SSHConnectService;
import org.forome.annotation.utils.ArgumentParser;
import org.forome.astorage.core.liftover.LiftoverConnector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//
//vulitin@ip-172-31-24-96:~$ PYTHONPATH=/data/bgm/versions/master/anfisa python -m annotations.singleton -a gnomad 1:103471457 "CCATCAT>CCAT"
//		Namespace(annotations='gnomad', input=['1:103471457', 'CCATCAT>CCAT'], test=1)
//		{
//		"exomes": {
//		"AC": 154564,
//		"AF": 0.8046520344841948,
//		"AN": 192088
//		},
//		"genomes": {
//		"AC": 23164,
//		"AF": 0.7641353829913571,
//		"AN": 30314
//		},
//		"overall": {
//		"AC": 177728,
//		"AF": 0.7991295042310771,
//		"AN": 222402
//		},
//		"popmax": "EAS",
//		"popmax_af": 0.9925577156743621,
//		"popmax_an": 13168,
//		"url": [
//		"http://gnomad.broadinstitute.org/variant/1-103471456-CCAT-C"
//		]
//		}
//		0.00427389144897
//


public class Service {

	private final static Logger log = LoggerFactory.getLogger(Service.class);

	private static Service instance = null;

	private final Thread.UncaughtExceptionHandler uncaughtExceptionHandler;
	private final QueryPool queryPool;

	private final Config config;
	private final ServiceConfig serviceConfig;
	private final SSHConnectService sshTunnelService;
	private final DatabaseConnectService databaseConnectService;
	private final DatabaseService databaseService;
	private final NetworkService networkService;

	private final SourceService sourceService;
	private final GnomadConnectorImpl gnomadConnector;
	private final SpliceAIConnector spliceAIConnector;
	private final HgmdConnector hgmdConnector;
	private final ClinvarConnector clinvarConnector;
	private final LiftoverConnector liftoverConnector;
//	private final FastaSourcePython fastaSource;
	private final GTFConnector gtfConnector;
	private final GTEXConnector gtexConnector;
	private final PharmGKBConnector pharmGKBConnector;
//	private final AStorageHttp sourceHttp38;
	private final EnsemblVepService ensemblVepService;
	private final AnfisaConnector anfisaConnector;


	private final NotificationService notificationService;

	public Service(ArgumentParser arguments, Thread.UncaughtExceptionHandler uncaughtExceptionHandler) throws Exception {
		instance = this;

		this.uncaughtExceptionHandler = uncaughtExceptionHandler;
		this.queryPool = new QueryPool(new QueryPoolExceptionBuilder() , uncaughtExceptionHandler);

		this.config = new Config();
		this.serviceConfig = new ServiceConfig();
		this.sshTunnelService = new SSHConnectService();
		this.databaseConnectService = new DatabaseConnectService(sshTunnelService, serviceConfig.databaseConfig);
		this.databaseService = new DatabaseService(this);
		this.networkService = new NetworkService(arguments.port, uncaughtExceptionHandler);

		if (serviceConfig.notificationSlackConfig != null) {
			this.notificationService = new NotificationService(serviceConfig.notificationSlackConfig);
		} else {
			this.notificationService = null;
		}

		this.sourceService = new SourceService(serviceConfig.sourceConfig);

//        this.gnomadConnector = new GnomadConnectorOld(databaseConnectService, serviceConfig.gnomadConfigConnector, uncaughtExceptionHandler);

		this.liftoverConnector = new LiftoverConnector();
//		this.fastaSource = new FastaSourcePython(databaseConnectService, serviceConfig.aStorageConfigConnector);

		this.gnomadConnector = new GnomadConnectorImpl(
				new GnomadDataSourceHttp(liftoverConnector, sourceService.dataSource),
				uncaughtExceptionHandler
		);
//		this.gnomadConnector = new GnomadConnectorImpl(databaseConnectService, serviceConfig.gnomadConfigConnector, uncaughtExceptionHandler);

		this.spliceAIConnector = new SpliceAIConnectorImpl(
				new SpliceAIDataSourceHttp(liftoverConnector)
		);
//		this.spliceAIConnector = new SpliceAIConnector(databaseConnectService, serviceConfig.spliceAIConfigConnector);

//		this.hgmdConnector = new HgmdConnectorHttp();
		this.hgmdConnector = new HgmdConnectorMysql(databaseConnectService, liftoverConnector, serviceConfig.hgmdConfigConnector);

//		this.clinvarConnector = new ClinvarConnectorHttp();
		this.clinvarConnector = new ClinvarConnectorMysql(databaseConnectService, liftoverConnector, serviceConfig.foromeConfigConnector);

//		this.gtfConnector = new GTFConnectorImpl(
//				new GTFDataSourceHttp(databaseConnectService, liftoverConnector, serviceConfig.aStorageConfigConnector),
//				uncaughtExceptionHandler
//		);
		this.gtfConnector = new GTFConnectorImpl(
				new GTFDataConnector(
						new DatabaseConnector(databaseConnectService, serviceConfig.gtfConfigConnector)
				),
				liftoverConnector,
				uncaughtExceptionHandler
		);

//		this.gtexConnector = new GTEXConnectorHttp();
		this.gtexConnector = new GTEXConnectorMysql(databaseConnectService, serviceConfig.foromeConfigConnector);

//		this.pharmGKBConnector = new PharmGKBConnectorHttp();
		this.pharmGKBConnector = new PharmGKBConnectorMysql(databaseConnectService, serviceConfig.foromeConfigConnector);

//		this.sourceHttp38 = new AStorageHttp(
//				databaseConnectService, liftoverConnector
//		);

		this.ensemblVepService = new EnsemblVepExternalService(uncaughtExceptionHandler);
//        this.ensemblVepService = new EnsemblVepInlineService(
//                sshTunnelService,
//                serviceConfig.ensemblVepConfigConnector,
//                new RefConnector(databaseConnectService, serviceConfig.refConfigConnector)
//        );

		this.anfisaConnector = new AnfisaConnector(
				sourceService,
				gnomadConnector,
				spliceAIConnector,
				hgmdConnector,
				clinvarConnector,
				liftoverConnector,
				gtfConnector,
				gtexConnector,
				pharmGKBConnector
//				,
//				sourceHttp38
		);

		queryPool.execute(this.databaseService.getDomainObjectSource(), new Query<Void>() {

			private ReadableResource<UserReadable> userReadableResource;
			private UserEditableComponent userEditableComponent;

			@Override
			public void prepare(ResourceProvider resources) {
				userReadableResource = resources.getReadableResource(UserReadable.class);
				userEditableComponent = new UserEditableComponent(resources);
			}

			@Override
			public Void execute(QueryTransaction transaction) throws AnnotatorException {
				if (!userReadableResource.iterator(transaction).hasNext()) {
					userEditableComponent.create("admin", "b82nfGl5sdg", transaction);
				}
				return null;
			}
		}).get();
	}

	public QueryPool getQueryPool() {
		return queryPool;
	}

	public Thread.UncaughtExceptionHandler getUncaughtExceptionHandler() {
		return uncaughtExceptionHandler;
	}

	public Config getConfig() {
		return config;
	}

	public ServiceConfig getServiceConfig() {
		return serviceConfig;
	}

	public DatabaseService getDatabaseService() {
		return databaseService;
	}

	public NetworkService getNetworkService() {
		return networkService;
	}

	public SourceService getSourceService() {
		return sourceService;
	}

	public DatabaseConnectService getDatabaseConnectService() {
		return databaseConnectService;
	}

	public GnomadConnectorImpl getGnomadConnector() {
		return gnomadConnector;
	}

	public LiftoverConnector getLiftoverConnector() {
		return liftoverConnector;
	}

	public GTFConnector getGtfConnector() {
		return gtfConnector;
	}

	public EnsemblVepService getEnsemblVepService() {
		return ensemblVepService;
	}

	public AnfisaConnector getAnfisaConnector() {
		return anfisaConnector;
	}

	public NotificationService getNotificationService() {
		return notificationService;
	}

	public void stop() {
		anfisaConnector.close();
		gtfConnector.close();
		liftoverConnector.close();
		clinvarConnector.close();
		hgmdConnector.close();
		spliceAIConnector.close();
		gnomadConnector.close();
		ensemblVepService.close();

		databaseConnectService.close();
		sshTunnelService.close();
		instance = null;
	}

	public static Service getInstance() {
		return instance;
	}
}
