/*
 *  Copyright (c) 2020. Vladimir Ulitin, Partners Healthcare and members of Forome Association
 *
 *  Developed by Vladimir Ulitin and Michael Bouzinier
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 * 	 http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package org.forome.annotation.makefavordatabase.validation.main;

import org.forome.database.exception.DatabaseException;
import org.forome.database.utils.TypeConvert;
import org.forome.rocksdb.RocksDBProvider;
import org.rocksdb.*;
import org.rocksdb.util.SizeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class ReadOnlyRocksDBConnector implements AutoCloseable {

	private final static Logger log = LoggerFactory.getLogger(ReadOnlyRocksDBConnector.class);

	public final Path pathRocksDB;

	public final RocksDB rocksDB;
	private final ConcurrentMap<String, ColumnFamilyHandle> columnFamilies;

	public ReadOnlyRocksDBConnector(Path pathRocksDB) throws Exception {
		this.pathRocksDB = pathRocksDB;

		log.debug("Load database: {}... ", pathRocksDB.toString());
		try (DBOptions options = buildOptions()) {
			List<ColumnFamilyDescriptor> columnFamilyDescriptors = getColumnFamilyDescriptors();

			List<ColumnFamilyHandle> columnFamilyHandles = new ArrayList<>();
			rocksDB = RocksDB.openReadOnly(options, pathRocksDB.toString(), columnFamilyDescriptors, columnFamilyHandles);

			columnFamilies = new ConcurrentHashMap<>();
			for (int i = 0; i < columnFamilyDescriptors.size(); i++) {
				String columnFamilyName = TypeConvert.unpackString(columnFamilyDescriptors.get(i).getName());
				ColumnFamilyHandle columnFamilyHandle = columnFamilyHandles.get(i);
				columnFamilies.put(columnFamilyName, columnFamilyHandle);
			}
		} catch (RocksDBException e) {
			log.error("Exception load database: {}", pathRocksDB.toString(), e);
			throw new DatabaseException(e);
		}
		log.debug("Load database: {}... complete", pathRocksDB.toString());
	}

	public ColumnFamilyHandle getColumnFamily(String name) {
		return columnFamilies.get(name);
	}

	public ColumnFamilyHandle createColumnFamily(String name) throws RocksDBException {
		ColumnFamilyDescriptor columnFamilyDescriptor = new ColumnFamilyDescriptor(name.getBytes(StandardCharsets.UTF_8));
		ColumnFamilyHandle columnFamilyHandle = rocksDB.createColumnFamily(columnFamilyDescriptor);
		columnFamilies.put(name, columnFamilyHandle);
		return columnFamilyHandle;
	}

	public void dropColumnFamily(String name) throws RocksDBException {
		rocksDB.dropColumnFamily(getColumnFamily(name));
		columnFamilies.remove(name);
	}

	private DBOptions buildOptions() throws RocksDBException {
		final String optionsFilePath = pathRocksDB.toString() + ".ini";

		DBOptions options = new DBOptions();
		if (Files.exists(Paths.get(optionsFilePath))) {
			final List<ColumnFamilyDescriptor> ignoreDescs = new ArrayList<>();
			OptionsUtil.loadOptionsFromFile(optionsFilePath, Env.getDefault(), options, ignoreDescs, false);
		} else {
			options
					.setInfoLogLevel(InfoLogLevel.WARN_LEVEL)
					.setMaxTotalWalSize(100L * SizeUnit.MB);
		}

		return options.setCreateIfMissing(true);
	}

	private List<ColumnFamilyDescriptor> getColumnFamilyDescriptors() throws RocksDBException {
		List<ColumnFamilyDescriptor> columnFamilyDescriptors = new ArrayList<>();

		try (Options options = new Options()) {
			for (byte[] columnFamilyName : RocksDB.listColumnFamilies(options, pathRocksDB.toString())) {
				columnFamilyDescriptors.add(new ColumnFamilyDescriptor(columnFamilyName));
			}
		}

		if (columnFamilyDescriptors.isEmpty()) {
			columnFamilyDescriptors.add(new ColumnFamilyDescriptor(TypeConvert.pack(RocksDBProvider.DEFAULT_COLUMN_FAMILY)));
		}

		return columnFamilyDescriptors;
	}

	@Override
	public void close() {
		rocksDB.close();
	}
}
