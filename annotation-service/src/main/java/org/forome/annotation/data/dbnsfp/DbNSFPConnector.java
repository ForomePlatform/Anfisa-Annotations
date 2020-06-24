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

package org.forome.annotation.data.dbnsfp;

import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import org.forome.annotation.data.anfisa.struct.AnfisaExecuteContext;
import org.forome.annotation.data.dbnsfp.struct.DbNSFPItem;
import org.forome.annotation.struct.variant.Variant;
import org.forome.annotation.utils.MathUtils;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class DbNSFPConnector {

	public List<DbNSFPItem> getAll(AnfisaExecuteContext context, Variant variant){
		JSONArray jRecords = (JSONArray) context.sourceSpliceAI_and_dbNSFP.get("dbNSFP");
		if (jRecords == null) {
			return Collections.emptyList();
		}

		List<JSONObject> records = jRecords.stream()
				.map(o -> (JSONObject) o)
				.filter(item -> item.getAsString("REF").equals(variant.getRef()) && item.getAsString("ALT").equals(variant.getStrAlt()))
				.collect(Collectors.toList());

		return records.stream().map(jsonObject -> _build(jsonObject)).collect(Collectors.toList());
	}

	private static DbNSFPItem _build(JSONObject jsonObject) {
		return new DbNSFPItem(
				MathUtils.toDouble(jsonObject.getAsNumber("CADD_phred"))
		);
	}
}
