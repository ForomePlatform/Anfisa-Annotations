package org.forome.annotation.connector.spliceai;

import org.forome.annotation.config.ServiceConfig;
import org.junit.Assert;
import org.junit.Before;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SpliceAIBaseTest {

	private final static Logger log = LoggerFactory.getLogger(SpliceAIBaseTest.class);

	protected SpliceAIConnector spliceAIConnector;

	@Before
	public void init() throws Throwable {
		ServiceConfig serviceConfig = new ServiceConfig();
		spliceAIConnector = new SpliceAIConnector(serviceConfig.spliceAIConfigConnector, (t, e) -> {
			log.error("Fail", e);
			Assert.fail();
		});
	}
}