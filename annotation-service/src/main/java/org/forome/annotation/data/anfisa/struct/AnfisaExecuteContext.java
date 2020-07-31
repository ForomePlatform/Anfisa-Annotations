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

package org.forome.annotation.data.anfisa.struct;

import net.minidev.json.JSONObject;
import org.forome.annotation.data.anfisa.AnfisaConnector;
import org.forome.annotation.data.astorage.struct.AStorageSource;
import org.forome.annotation.data.fasta.FastaSource;
import org.forome.annotation.struct.Assembly;
import org.forome.annotation.struct.Interval;
import org.forome.annotation.struct.Sequence;
import org.forome.annotation.struct.variant.Variant;

import java.util.HashMap;
import java.util.Map;

public class AnfisaExecuteContext {

	private final String CACHE_MASKED_REGION = "masked_region";

	public final AnfisaInput anfisaInput;

	public final Variant variant;
	public final JSONObject vepJson;

	public Double gnomadAfFam;

	public AStorageSource sourceAStorageHttp;

	private final Map<String, Object> cache;

	public AnfisaExecuteContext(
			AnfisaInput anfisaInput,
			Variant variant,
			JSONObject vepJson
	) {
		this.anfisaInput = anfisaInput;

		this.variant = variant;
		this.vepJson = vepJson;

		this.cache = new HashMap<>();
	}

	public boolean getMaskedRegion(AnfisaConnector anfisaConnector, AnfisaExecuteContext context) {
		return (boolean) cache.computeIfAbsent(CACHE_MASKED_REGION, s -> {
			FastaSource fastaSource = anfisaConnector.fastaSource;
			Assembly assembly = anfisaInput.mCase.assembly;
			Interval interval = Interval.of(
					variant.chromosome,
					variant.getStart(),
					(variant.getStart() < variant.end) ? variant.end : variant.getStart()
			);
			Sequence sequence = fastaSource.getSequence(context, assembly, interval);
			String vSequence = sequence.value;

			//Если есть маленькие буквы, то мы имеем дело с замаскированными регионами тандемных повторов
			boolean maskedRegion = !vSequence.equals(vSequence.toUpperCase());

			return maskedRegion;
		});
	}
}
