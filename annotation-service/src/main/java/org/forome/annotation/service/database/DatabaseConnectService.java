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

package org.forome.annotation.service.database;

import com.infomaximum.database.exception.DatabaseException;
import com.mchange.v2.c3p0.ComboPooledDataSource;
import org.forome.annotation.config.connector.base.DatabaseConfigConnector;
import org.forome.annotation.config.database.DatabaseConfig;
import org.forome.annotation.config.sshtunnel.SshTunnelConfig;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.service.database.rocksdb.annotator.SourceDatabase;
import org.forome.annotation.service.database.rocksdb.favor.FavorDatabase;
import org.forome.annotation.service.ssh.SSHConnectService;
import org.forome.annotation.service.ssh.struct.SSHConnect;
import org.forome.annotation.struct.Assembly;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.time.Duration;
import java.util.HashMap;
import java.util.Map;

public class DatabaseConnectService implements AutoCloseable {

	private final static Logger log = LoggerFactory.getLogger(DatabaseConnectService.class);

	private final SourceDatabase sourceDatabase37;
	private final SourceDatabase sourceDatabase38;

	private final FavorDatabase favorDatabase;

	private final SSHConnectService sshTunnelService;
	private final Map<String, ComboPooledDataSource> dataSources;

	public DatabaseConnectService(SSHConnectService sshTunnelService, DatabaseConfig databaseConfig) throws DatabaseException {
		this.sshTunnelService = sshTunnelService;
		this.dataSources = new HashMap<>();

		if (databaseConfig.hg37 != null) {
			sourceDatabase37 = new SourceDatabase(Assembly.GRCh37, databaseConfig.hg37);
		} else {
			sourceDatabase37 = null;
		}
		if (databaseConfig.hg38 != null) {
			sourceDatabase38 = new SourceDatabase(Assembly.GRCh38, databaseConfig.hg38);
		} else {
			sourceDatabase38 = null;
		}
		if (databaseConfig.favor != null) {
			favorDatabase = new FavorDatabase(databaseConfig.favor);
		} else {
			favorDatabase = null;
		}
	}

	public ComboPooledDataSource getDataSource(DatabaseConfigConnector databaseConfigConnector) {
		String keyDataSource = getKeyDataSource(databaseConfigConnector);
		ComboPooledDataSource dataSource = dataSources.get(keyDataSource);
		if (dataSource == null) {
			synchronized (dataSources) {
				dataSource = dataSources.get(keyDataSource);
				if (dataSource == null) {
					dataSource = connect(databaseConfigConnector);
					dataSources.put(keyDataSource, dataSource);
				}
			}
		}
		return dataSource;
	}

	private ComboPooledDataSource connect(DatabaseConfigConnector databaseConfigConnector) {
		try {
			StringBuilder jdbcUrl = new StringBuilder("jdbc:mysql://")
					.append(databaseConfigConnector.mysqlHost).append(':');

			int mysqlUrlPort;
			SshTunnelConfig sshTunnelConfig = databaseConfigConnector.sshTunnelConfig;
			if (sshTunnelConfig != null) {
				SSHConnect sshTunnel = sshTunnelService.getSSHConnect(
						sshTunnelConfig.host,
						sshTunnelConfig.port,
						sshTunnelConfig.user,
						sshTunnelConfig.key
				);
				mysqlUrlPort = sshTunnel.getTunnel(databaseConfigConnector.mysqlPort);
			} else {
				mysqlUrlPort = databaseConfigConnector.mysqlPort;
			}
			jdbcUrl.append(mysqlUrlPort).append('/').append(databaseConfigConnector.mysqlDatabase)
					.append("?user=").append(databaseConfigConnector.mysqlUser)
					.append("&password=").append(databaseConfigConnector.mysqlPassword)
					.append("&useSSL=false");

			String driverName = "com.mysql.jdbc.Driver";
			Class.forName(driverName).newInstance();

			ComboPooledDataSource pooledDataSource = new ComboPooledDataSource();
			pooledDataSource.setDriverClass(driverName);
			pooledDataSource.setJdbcUrl(jdbcUrl.toString());
			pooledDataSource.setMinPoolSize(1);
			pooledDataSource.setAcquireIncrement(1);
			pooledDataSource.setMaxPoolSize(20);
			pooledDataSource.setCheckoutTimeout((int) Duration.ofMinutes(1).toMillis());
			pooledDataSource.setTestConnectionOnCheckin(false);
			pooledDataSource.setTestConnectionOnCheckout(true);

			log.debug("Database connected by: {}", databaseConfigConnector.mysqlHost);

			return pooledDataSource;
		} catch (Throwable ex) {
			throw ExceptionBuilder.buildExternalDatabaseException(ex);
		}
	}

	private static String getKeyDataSource(DatabaseConfigConnector databaseConfigConnector) {
		StringBuilder builderKey = new StringBuilder();
		if (databaseConfigConnector.sshTunnelConfig != null) {
			SshTunnelConfig sshTunnelConfig = databaseConfigConnector.sshTunnelConfig;
			String keySSHConnect = SSHConnectService.getKeySSHConnect(
					sshTunnelConfig.host, sshTunnelConfig.port, sshTunnelConfig.user
			);
			builderKey.append(keySSHConnect);
		}
		return builderKey.append(databaseConfigConnector.mysqlHost)
				.append(databaseConfigConnector.mysqlPort)
				.append(databaseConfigConnector.mysqlUser)
				.toString();
	}

	@Override
	public void close() {
		for (ComboPooledDataSource dataSource : dataSources.values()) {
			dataSource.close();
		}
	}

	public Source getSource(Assembly assembly) {
		switch (assembly) {
			case GRCh37:
				return sourceDatabase37;
			case GRCh38:
				return sourceDatabase38;
			default:
				throw new RuntimeException();
		}
	}

	public FavorDatabase getFavorDatabase() {
		return favorDatabase;
	}
}
