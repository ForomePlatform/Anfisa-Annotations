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

package org.forome.annotation.processing.smavariant;

import org.forome.annotation.processing.smavariant.vcf.SplitMAVariantVcf;
import org.forome.annotation.struct.mavariant.MAVariant;
import org.forome.annotation.struct.mavariant.MAVariantVCF;
import org.forome.annotation.struct.variant.Variant;

import java.util.List;

public abstract class SplitMAVariant {

	public abstract List<Variant> split();

	public static SplitMAVariant build(MAVariant maVariant) {
		if (maVariant instanceof MAVariantVCF) {
			return new SplitMAVariantVcf((MAVariantVCF) maVariant);
		} else {
			throw new RuntimeException();
		}
	}
}

