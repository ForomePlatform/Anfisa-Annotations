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

package org.forome.annotation.service.source.internal;

import org.forome.database.exception.DatabaseException;
import org.forome.annotation.config.source.SourceInternalConfig;
import org.forome.annotation.service.source.DataSource;
import org.forome.annotation.service.source.external.ExternalDataSource;
import org.forome.annotation.service.source.external.source.ExternalSource;
import org.forome.annotation.service.source.internal.source.InternalSource;
import org.forome.annotation.service.source.struct.Source;
import org.forome.astorage.AStorage;
import org.forome.core.struct.Assembly;

public class InternalDataSource implements DataSource {

	private final AStorage aStorage;

	private final ExternalDataSource httpDataSource;

	public InternalDataSource(SourceInternalConfig config) throws DatabaseException {
		AStorage.Builder builder = new AStorage.Builder();
		if (config.hg37 != null) {
			builder.withSource(Assembly.GRCh37, config.hg37);
		}
		if (config.hg38 != null) {
			builder.withSource(Assembly.GRCh38, config.hg38);
		}
		if (config.pastorage != null) {
			builder.withSourcePAStorage(config.pastorage);
		}
		aStorage = builder.build();

		httpDataSource = new ExternalDataSource(config.sourceExternalConfig);
	}

	public AStorage getAStorage() {
		return aStorage;
	}

	@Override
	public Source getSource(Assembly assembly) {
		switch (assembly) {
			case GRCh37:
				return new InternalSource(assembly, aStorage.getPaStorage(), aStorage.sourceDatabase37, (ExternalSource)httpDataSource.getSource(assembly));
			case GRCh38:
				return new InternalSource(assembly, aStorage.getPaStorage(), aStorage.sourceDatabase38, (ExternalSource)httpDataSource.getSource(assembly));
			default:
				throw new RuntimeException();
		}
	}


}
