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

package org.forome.annotation.annotator.utils;

import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import net.minidev.json.parser.ParseException;
import org.forome.annotation.struct.mcase.Cohort;
import org.forome.annotation.struct.mcase.MCase;
import org.forome.annotation.struct.mcase.Sample;
import org.forome.annotation.struct.mcase.Sex;
import org.forome.core.struct.Assembly;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class CaseUtils {

	/**
	 * sample["name"]      = sample_map.get(id, id)
	 * sample['family']    = family
	 * sample['id']        = id
	 * sample['father']    = father
	 * sample['mother']    = mother
	 * sample['sex']       = int(sex)
	 * sample['affected']  = (int(affected) == 2)
	 */
	public static MCase parseFamFile(Assembly assembly, InputStream isFam, InputStream isPatientIds, InputStream isCohorts) throws IOException, ParseException {
		LinkedHashMap<String, Sample> samples = new LinkedHashMap<>();

		Map<String, String> sampleNameMap = new HashMap<>();
		if (isPatientIds != null) {
			try (BufferedReader isBFamSampleName = new BufferedReader(new InputStreamReader(isPatientIds))) {
				String line;
				while ((line = isBFamSampleName.readLine()) != null) {
					line = line.trim();
					if (line.length() == 0) continue;

					String[] values = Arrays.stream(line.split("\t")).filter(s -> !s.trim().isEmpty()).toArray(String[]::new);
					if (values.length != 4) {
						throw new RuntimeException("Not support FamSampleName format");
					}
					String id = values[3];
					String name = values[0];
					sampleNameMap.put(id, name);
				}
			}
		}

		//Раскидываем маски семплов по когортам
		LinkedHashMap<Cohort, List<Pattern>> cohortMaskSamples = new LinkedHashMap<>();
		if (isCohorts != null) {
			JSONArray jCohorts = (JSONArray) new JSONParser(JSONParser.DEFAULT_PERMISSIVE_MODE).parse(isCohorts);
			for (Object joCohort : jCohorts) {
				JSONObject jCohort = (JSONObject) joCohort;
				Cohort cohort = new Cohort(jCohort.getAsString("name"));
				List<Pattern> maskSamples = ((JSONArray) jCohort.get("members")).stream().map(o -> (String) o)
						.map(s -> Pattern.compile(
								"^" + s.replaceAll("\\*", "(.*)") + "$"
						))
						.collect(Collectors.toList());
				cohortMaskSamples.put(cohort, maskSamples);
			}
			isCohorts.close();
		}

		try (BufferedReader isBFam = new BufferedReader(new InputStreamReader(isFam))) {
			String line;
			while ((line = isBFam.readLine()) != null) {
				if (line.startsWith("#")) {
					continue;
				}
				String[] sl = Arrays.stream(
						line.split("\\s+")
				).filter(s -> !s.isEmpty()).toArray(String[]::new);

				String id = sl[1];
				String name = sampleNameMap.getOrDefault(id, id);
				String family = sl[0];
				String father = sl[2];
				String mother = sl[3];
				int isex = Integer.parseInt(sl[4]);
				boolean affected = Integer.parseInt(sl[5]) == 2;

				Cohort cohort = null;
				for (Map.Entry<Cohort, List<Pattern>> entry : cohortMaskSamples.entrySet()) {
					for (Pattern pattern : entry.getValue()) {
						if (pattern.matcher(id).matches()) {
							if (cohort != null) {
								throw new RuntimeException("Cohort: Not unique matches sample: " + id);
							}

							cohort = entry.getKey();
							break;
						}
					}
				}

				Sex sex;
				if (isex == 0 ) {
					sex = Sex.UNKNOWN;
				} else if (isex == 1 ) {
					sex = Sex.MALE;
				} else if (isex == 2 ) {
					sex = Sex.FEMALE;
				} else {
					throw new RuntimeException("Unknown sex: " + isex);
				}

				Sample sample = new Sample(id, name, family, father, mother, sex, affected, cohort);
				samples.put(id, sample);
			}
		}

		return new MCase.Builder(assembly, samples, new ArrayList<>(cohortMaskSamples.keySet())).build();
	}

}
