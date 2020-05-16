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

package org.forome.annotation.data.gtex.http;

import org.forome.annotation.data.gtex.GTEXConnector;
import org.forome.annotation.data.gtex.struct.Tissue;
import org.forome.annotation.struct.SourceMetadata;

import java.util.Collections;
import java.util.List;

public class GTEXConnectorHttp implements GTEXConnector {

	@Override
	public List<SourceMetadata> getSourceMetadata() {
		return Collections.emptyList();
	}

	@Override
	public List<Tissue> getTissues(String gene) {
		return Collections.emptyList();
	}
}
