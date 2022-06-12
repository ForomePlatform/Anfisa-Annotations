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

package org.forome.annotation.service.source;

import org.forome.database.exception.DatabaseException;
import org.forome.annotation.config.source.SourceConfig;
import org.forome.annotation.service.source.external.ExternalDataSource;
import org.forome.annotation.service.source.internal.InternalDataSource;
import org.forome.annotation.service.source.wrapper.WrapperDataSource;

public class SourceService {

	public final WrapperDataSource dataSource;

	public SourceService(SourceConfig config) {

		DataSource ds;
		if (config.sourceInternalConfig != null) {
			try {
				ds = new InternalDataSource(config.sourceInternalConfig);
			} catch (DatabaseException e) {
				throw new RuntimeException(e);
			}
		} else if (config.sourceExternalConfig != null) {
			ds = new ExternalDataSource(config.sourceExternalConfig);
		} else {
			throw new RuntimeException();
		}

		dataSource = new WrapperDataSource(ds);
	}

}
