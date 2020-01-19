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

package org.forome.annotation.processing.graphql.record.view.bioinformatics;

import graphql.annotations.annotationTypes.GraphQLField;
import graphql.annotations.annotationTypes.GraphQLName;
import org.forome.annotation.struct.Allele;
import org.forome.annotation.struct.mcase.MCase;
import org.forome.annotation.struct.mcase.Sex;
import org.forome.annotation.struct.variant.Variant;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

@GraphQLName("record_view_bioinformatics")
public class GRecordViewBioinformatics {

	public final MCase mCase;
	public final Variant variant;

	public GRecordViewBioinformatics(MCase mCase, Variant variant) {
		this.variant = variant;
		this.mCase = mCase;
	}

	@GraphQLField
	@GraphQLName("zygosity")
	public String getZygosity() {
		if (mCase == null) {
			return null;
		}

		List<Allele> alleles = variant.getGenotype(mCase.proband).getAllele();
		if (alleles == null) {
			return null;
		}

		Set<Allele> uniqueAllelies = new HashSet<>(alleles);

		/**
		 * Ситуация когда ни один аллель генотипа не относится к иследуемому варианту,
		 * например ситуация, когда мы разрезали мультиалельный вариант
		 */
		if (!uniqueAllelies.contains(variant.getRefAllele()) && !uniqueAllelies.contains(variant.getAlt())) {
			return "Unknown";
		}

		String chr = variant.chromosome.getChar();
		if ("X".equals(chr.toUpperCase()) && mCase.proband.sex == Sex.MALE) {
			return "X-linked";
		}
		if (uniqueAllelies.size() == 1) {
			return "Homozygous";
		}
		if (uniqueAllelies.size() == 2) {
			return "Heterozygous";
		}
		return "Unknown";
	}
}
