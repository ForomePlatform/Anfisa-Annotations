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

package org.forome.annotation.service.database.rocksdb.favor;

import org.forome.database.exception.DatabaseException;
import org.forome.annotation.utils.compression.GZIPCompression;
import org.forome.astorage.core.rocksdb.RocksDBDatabase;
import org.forome.astorage.core.utils.bits.IntegerBits;
import org.rocksdb.ColumnFamilyHandle;
import org.rocksdb.RocksDBException;
import org.rocksdb.RocksIterator;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FavorDatabase extends RocksDBDatabase {

	public static final String COLUMN_FAMILY_DATA = "data";
//	public static final String COLUMN_FAMILY_META = "meta";

	public FavorDatabase(Path pathDatabase) throws DatabaseException {
		super(pathDatabase);
	}

	public int getSize() {
		try (RocksIterator rocksIterator = rocksDB.newIterator(getDataColumnFamily())) {
			rocksIterator.seekToLast();

			byte[] bytes = rocksIterator.key();
			return getOrder(bytes) + 1;
		}
	}

	public String getRecord(int order) throws RocksDBException {
		byte[] key = getKeyData(order);
		byte[] value = rocksDB.get(getDataColumnFamily(), key);
		if (value == null) return null;

		return GZIPCompression.decompress(value);
	}

	public List<String> getSequenceRecords(int minOrder, int maxOrder) {
		List<String> records = new ArrayList<>();
		try (RocksIterator rocksIterator = rocksDB.newIterator(getDataColumnFamily())) {
			rocksIterator.seek(getKeyData(minOrder));
			while (rocksIterator.isValid()) {
				int order = FavorDatabase.getOrder(rocksIterator.key());
				if (order > maxOrder) break;

				String record = GZIPCompression.decompress(rocksIterator.value());
				records.add(record);

				rocksIterator.next();
			}
		}

		return records;
	}

	private ColumnFamilyHandle getDataColumnFamily() {
		return getColumnFamily(COLUMN_FAMILY_DATA);
	}

//	private ColumnFamilyHandle getMetaColumnFamily() {
//		return getColumnFamily(COLUMN_FAMILY_META);
//	}

	public static byte[] getKeyData(int order) {
		return IntegerBits.toByteArray(order);
	}

	public static int getOrder(byte[] key) {
		return IntegerBits.fromByteArray(key);
	}

//	public static byte[] getKeyMetaSize() {
//		return StringBits.toByteArray("size");
//	}
}
