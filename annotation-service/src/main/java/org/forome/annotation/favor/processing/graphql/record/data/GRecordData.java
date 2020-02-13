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

package org.forome.annotation.favor.processing.graphql.record.data;

import graphql.annotations.annotationTypes.GraphQLField;
import graphql.annotations.annotationTypes.GraphQLName;
import org.forome.annotation.favor.utils.struct.table.Row;

@GraphQLName("record_data")
public class GRecordData {

	public final Row row;

	public GRecordData(Row row) {
		this.row = row;
	}

	@GraphQLField
	@GraphQLName("label")
	public String getLabel() {
		return row.getValue("variant_format");
	}

	@GraphQLField
	@GraphQLName("color_code")
	public String getColorCode() {
		return null;
	}
}
