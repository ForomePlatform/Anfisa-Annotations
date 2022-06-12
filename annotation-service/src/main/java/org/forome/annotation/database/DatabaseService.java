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

package org.forome.annotation.database;

import org.forome.database.domainobject.DomainObjectSource;
import org.forome.database.maintenance.ChangeMode;
import org.forome.database.maintenance.SchemaService;
import org.forome.database.provider.DBProvider;
import org.forome.database.schema.Schema;
import org.forome.rocksdb.RocksDataBaseBuilder;
import org.forome.annotation.Service;
import org.forome.annotation.database.entityobject.user.UserReadable;

import java.nio.file.Path;

public class DatabaseService {

	private final Service service;

	private DomainObjectSource domainObjectSource;

	public DatabaseService(Service service) throws Exception {
		this.service = service;

		Path dataPath = service.getConfig().dataPath;

		DBProvider dbProvider = new RocksDataBaseBuilder()
				.withPath(dataPath.resolve("database"))
				.build();

		Schema schema;
		if (Schema.exists(dbProvider)) {
			schema = Schema.read(dbProvider);
		} else {
			schema = Schema.create(dbProvider);
		}
		Schema.resolve(UserReadable.class);

		this.domainObjectSource = new DomainObjectSource(dbProvider);

		new SchemaService(dbProvider)
				.setNamespace("org.forome.annotation")
				.setChangeMode(ChangeMode.CREATION)
				.setSchema(schema)
				.execute();
	}

	public DomainObjectSource getDomainObjectSource() {
		return domainObjectSource;
	}
}
