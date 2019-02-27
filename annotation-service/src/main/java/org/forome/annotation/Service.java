package org.forome.annotation;

import org.forome.annotation.config.ServiceConfig;
import org.forome.annotation.connector.anfisa.AnfisaConnector;
import org.forome.annotation.connector.clinvar.ClinvarConnector;
import org.forome.annotation.connector.gnomad.GnomadConnector;
import org.forome.annotation.connector.hgmd.HgmdConnector;
import org.forome.annotation.connector.liftover.LiftoverConnector;
import org.forome.annotation.exception.ServiceException;
import org.forome.annotation.executionqueue.*;
import org.forome.annotation.network.NetworkService;
import org.forome.annotation.network.component.UserEditableComponent;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.forome.annotation.connector.gtf.GTFConnector;
import org.forome.annotation.database.DatabaseService;
import org.forome.annotation.database.entityobject.user.UserReadable;
import org.forome.annotation.utils.ArgumentParser;

import java.io.IOException;

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
	private final ExecutionQueue executionQueue;

	private final ServiceConfig serviceConfig;
	private final DatabaseService databaseService;
	private final NetworkService networkService;

	private final GnomadConnector gnomadConnector;
	private final HgmdConnector hgmdConnector;
	private final ClinvarConnector clinvarConnector;
	private final LiftoverConnector liftoverConnector;
	private final GTFConnector gtfConnector;
	private final AnfisaConnector anfisaConnector;

	public Service(ArgumentParser arguments, Thread.UncaughtExceptionHandler uncaughtExceptionHandler) throws Exception {
		instance = this;

		this.uncaughtExceptionHandler = uncaughtExceptionHandler;
		this.executionQueue = new ExecutionQueue(uncaughtExceptionHandler);

		this.serviceConfig = new ServiceConfig();
		this.databaseService = new DatabaseService(this);
		this.networkService = new NetworkService(arguments.port, uncaughtExceptionHandler);

		this.gnomadConnector = new GnomadConnector(serviceConfig.gnomadConfigConnector, uncaughtExceptionHandler);
		this.hgmdConnector = new HgmdConnector(serviceConfig.hgmdConfigConnector);
		this.clinvarConnector = new ClinvarConnector(serviceConfig.clinVarConfigConnector);
		this.liftoverConnector = new LiftoverConnector();
		this.gtfConnector = new GTFConnector(serviceConfig.gtfConfigConnector, uncaughtExceptionHandler);
		this.anfisaConnector = new AnfisaConnector(gnomadConnector, hgmdConnector, clinvarConnector, liftoverConnector);

		executionQueue.execute(this, new Execution<Void>() {

			private ReadableResource<UserReadable> userReadableResource;
			private UserEditableComponent userEditableComponent;

			@Override
			public void prepare(ResourceProvider resources) {
				userReadableResource = resources.getReadableResource(UserReadable.class);
				userEditableComponent = new UserEditableComponent(resources);
			}

			@Override
			public Void execute(ExecutionTransaction transaction) throws ServiceException {
				if (!userReadableResource.iterator(transaction).hasNext()) {
					userEditableComponent.create("admin", "b82nfGl5sdg", transaction);
				}
				return null;
			}
		}).get();
	}

	public ExecutionQueue getExecutionQueue() {
		return executionQueue;
	}

	public Thread.UncaughtExceptionHandler getUncaughtExceptionHandler() {
		return uncaughtExceptionHandler;
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

	public GnomadConnector getGnomadConnector() {
		return gnomadConnector;
	}

	public LiftoverConnector getLiftoverConnector() {
		return liftoverConnector;
	}

	public GTFConnector getGtfConnector() {
		return gtfConnector;
	}

	public AnfisaConnector getAnfisaConnector() {
		return anfisaConnector;
	}

	public void stop() {
		try {
			gnomadConnector.close();
			anfisaConnector.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		instance = null;
	}

	public static Service getInstance() {
		return instance;
	}
}
