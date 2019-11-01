package org.forome.annotation.connector.gtex;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import org.forome.annotation.config.connector.GTEXConfigConnector;
import org.forome.annotation.connector.DatabaseConnector;
import org.forome.annotation.connector.gtex.struct.Tissue;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.service.database.DatabaseConnectService;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.concurrent.ExecutionException;

public class GTEXConnector implements AutoCloseable {

	private final DatabaseConnector databaseConnector;

	private final Map<Integer, String> tissueTypes;

	private final Cache<String, List<Tissue>> cacheTissues;

	public GTEXConnector(
			DatabaseConnectService databaseConnectService,
			GTEXConfigConnector gtexConfigConnector
	) throws Exception {
		this.databaseConnector = new DatabaseConnector(databaseConnectService, gtexConfigConnector);

		this.tissueTypes = Collections.unmodifiableMap(buildTissueTypes(databaseConnector));

		this.cacheTissues = CacheBuilder.newBuilder()
				.maximumSize(1000)
				.build();
	}

	public List<DatabaseConnector.Metadata> getMetadata() {
		return databaseConnector.getMetadata();
	}

	public List<Tissue> getTissues(String gene) {
		try {
			return cacheTissues.get(gene, () -> Collections.unmodifiableList(loadTissues(gene)));
		} catch (ExecutionException e) {
			Throwable cause = e.getCause();
			if (cause instanceof RuntimeException) {
				throw (RuntimeException) cause;
			} else {
				throw new RuntimeException(cause);
			}
		}
	}

	private List<Tissue> loadTissues(String gene) {
		String sql = String.format(
				"select TissueNo, Expression, RelExp from %s.GENE2TISSUE where GeneName = (select GeneName from %s.GENE where Description='%s')",
				databaseConnector.getDatabase(),
				databaseConnector.getDatabase(),
				gene
		);
		List<Tissue> tissues = new ArrayList<>();
		try (Connection connection = databaseConnector.createConnection()) {
			try (Statement statement = connection.createStatement()) {
				try (ResultSet resultSet = statement.executeQuery(sql)) {
					while (resultSet.next()) {
						String tissueName = tissueTypes.get(resultSet.getInt("TissueNo"));
						Tissue tissue = new Tissue(
								tissueName,
								resultSet.getFloat("Expression"),
								resultSet.getFloat("RelExp")
						);
						tissues.add(tissue);
					}
				}
			}
		} catch (SQLException ex) {
			throw ExceptionBuilder.buildExternalDatabaseException(ex);
		}
		return tissues;
	}


	private static Map<Integer, String> buildTissueTypes(DatabaseConnector databaseConnector) {
		String sql = String.format(
				"SELECT TissueNo, Name FROM %s.TISSUE", databaseConnector.getDatabase()
		);
		Map<Integer, String> tissueTypes = new HashMap<>();
		try (Connection connection = databaseConnector.createConnection()) {
			try (Statement statement = connection.createStatement()) {
				try (ResultSet resultSet = statement.executeQuery(sql)) {
					while (resultSet.next()) {
						tissueTypes.put(
								resultSet.getInt("TissueNo"),
								resultSet.getString("Name")
						);
					}
				}
			}
		} catch (SQLException ex) {
			throw ExceptionBuilder.buildExternalDatabaseException(ex);
		}
		return tissueTypes;
	}

	@Override
	public void close() {
		this.databaseConnector.close();
	}
}
