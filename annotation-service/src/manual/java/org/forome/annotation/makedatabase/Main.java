/*
 *  Copyright (c) 2019. Vladimir Ulitin, Partners Healthcare and members of Forome Association
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

package org.forome.annotation.makedatabase;

import org.forome.database.utils.TypeConvert;
import org.forome.rocksdb.RocksDBProvider;
import org.forome.astorage.core.packer.PackInterval;
import org.forome.astorage.core.source.SourceDatabase;
import org.forome.core.struct.Interval;
import org.rocksdb.*;
import org.rocksdb.util.SizeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Main {

	private final static Logger log = LoggerFactory.getLogger(Main.class);

	private static RocksDB rocksDB;
	private static Map<String, ColumnFamilyHandle> columnFamilies;

	public static void main(String[] args) throws Exception {

		Path path = Paths.get("/home/kris/processtech/annotation-database");

		try (DBOptions options = buildOptions(path)) {
			List<ColumnFamilyDescriptor> columnFamilyDescriptors = getColumnFamilyDescriptors(path);

			List<ColumnFamilyHandle> columnFamilyHandles = new ArrayList<>();
			rocksDB = RocksDB.openReadOnly(options, path.toString(), columnFamilyDescriptors, columnFamilyHandles);

			columnFamilies = new HashMap<>();
			for (int i = 0; i < columnFamilyDescriptors.size(); i++) {
				String columnFamilyName = TypeConvert.unpackString(columnFamilyDescriptors.get(i).getName());
				ColumnFamilyHandle columnFamilyHandle = columnFamilyHandles.get(i);
				columnFamilies.put(columnFamilyName, columnFamilyHandle);
			}


			//
			int size = 0;
			int size4 = 0;
			int size16 = 0;
			int size100 = 0;
			try (RocksIterator rocksIterator = rocksDB.newIterator(
					columnFamilies.get(SourceDatabase.COLUMN_FAMILY_RECORD)
			)) {
				rocksIterator.seekToFirst();
				while (rocksIterator.isValid()) {

					Interval interval = new PackInterval().fromByteArray(rocksIterator.key());
					if (interval.start % 10000000 == 0) {
						log.debug("interval, chr: {}, pos: {}", interval.chromosome, interval.start);
					}

					ByteBuffer byteBuffer = ByteBuffer.wrap(rocksIterator.value());
					ArrayList<Short> numbers = new ArrayList<>();
					for (int i = 0; i < 100; i++) {
						numbers.add(byteBuffer.getShort());
						numbers.add(byteBuffer.getShort());
					}
					Set<Short> unique = new HashSet<>(numbers);

					size++;
					if (unique.size() < 4) {
						size4++;
					} else if (unique.size() < 16) {
						size16++;
					} else if (unique.size() < 100) {
						size100++;
					}
					rocksIterator.next();
				}
			}

			log.debug("size: {}", size);
			log.debug("<4, size4: {} vs {}% ", size4, size4 * 100 / size);
			log.debug("<16, size16: {} vs {}% ", size16, size16 * 100 / size);
			log.debug("<100, size100: {} vs {}% ", size100, size100 * 100 / size);

			System.exit(0);
		} catch (RocksDBException e) {
			log.debug("Exception: ", e);
			System.exit(1);
		}
	}


	private static DBOptions buildOptions(Path pathDatabase) throws RocksDBException {
		final String optionsFilePath = pathDatabase.toString() + ".ini";

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

	private static List<ColumnFamilyDescriptor> getColumnFamilyDescriptors(Path pathDatabase) throws RocksDBException {
		List<ColumnFamilyDescriptor> columnFamilyDescriptors = new ArrayList<>();

		try (Options options = new Options()) {
			for (byte[] columnFamilyName : RocksDB.listColumnFamilies(options, pathDatabase.toString())) {
				columnFamilyDescriptors.add(new ColumnFamilyDescriptor(columnFamilyName));
			}
		}

		if (columnFamilyDescriptors.isEmpty()) {
			columnFamilyDescriptors.add(new ColumnFamilyDescriptor(TypeConvert.pack(RocksDBProvider.DEFAULT_COLUMN_FAMILY)));
		}

		return columnFamilyDescriptors;
	}

}
