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

package org.forome.annotation.favor.processing.graphql.record.filters;

import graphql.annotations.annotationTypes.GraphQLField;
import graphql.annotations.annotationTypes.GraphQLName;
import org.forome.annotation.favor.utils.struct.table.Row;

import java.util.Arrays;

@GraphQLName("record_filters")
public class GRecordFilters {

	public final Row row;

	private final String[] splitvariantFormat;

	public GRecordFilters(Row row) {
		this.row = row;

		String variantFormat = row.getValue("variant_format");
		splitvariantFormat = variantFormat.split("-");
		if (splitvariantFormat.length != 4) throw new RuntimeException();
	}

	@GraphQLField
	@GraphQLName("type")
	public String getType() {
		String type = row.getValue("Type");
		return (type == null) ? "" : type;
	}

	@GraphQLField
	@GraphQLName("chromosome")
	public String getChromosome() {
		return splitvariantFormat[0];
	}

	@GraphQLField
	@GraphQLName("start")
	public int getStart() {
		return Integer.parseInt(splitvariantFormat[1]);
	}

	@GraphQLField
	@GraphQLName("end")
	public int getEnd() {
		return Integer.parseInt(splitvariantFormat[1]);
	}

	@GraphQLField
	@GraphQLName("ref")
	public String getRef() {
		return splitvariantFormat[2];
	}

	@GraphQLField
	@GraphQLName("alt")
	public String getAlt() {
		return splitvariantFormat[3];
	}

	@GraphQLField
	@GraphQLName("genes")
	public String[] getGenes() {
		return new String[]{ row.getValue("GeneName") };
	}

	@GraphQLField
	@GraphQLName("gnomad_total_af")
	public double getGnomadTotalAf() {
		String value = row.getValue("GNOMAD_total");
		if (value == null) return 0;
		return Double.parseDouble(value);
	}

	@GraphQLField
	@GraphQLName("gencode_category")
	public String getGencodeCategory() {
		return row.getValue("GENCODE.Category");
	}

	@GraphQLField
	@GraphQLName("gencode_exonic_category")
	public String getGencodeExonicCategory() {
		return row.getValue("GENCODE.EXONIC.Category");
	}

	@GraphQLField
	@GraphQLName("polyphen2_hdiv")
	public String getPolyphen2HDIV() {
		String value = row.getValue("Polyphen2_HDIV_pred");
		if (".".equals(value)) return null;
		return value;
	}

	@GraphQLField
	@GraphQLName("polyphen2_hvar")
	public String getPolyphen2HVAR() {
		String value = row.getValue("Polyphen2_HVAR_pred");
		if (".".equals(value)) return null;
		return value;
	}

	@GraphQLField
	@GraphQLName("polyphen_cat")
	public String getPolyPhenCat() {
		return row.getValue("PolyPhenCat");
	}

	@GraphQLField
	@GraphQLName("sift_cat")
	public String getSiftCat() {
		return row.getValue("SIFTcat");
	}

	@GraphQLField
	@GraphQLName("clinvar")
	public String[] getClinvar() {
		String value = row.getValue("clinvar");
		if (value==null) return null;
		return Arrays.stream(value.split("\\|"))
				.map(s -> s.trim())
				.filter(s -> !s.isEmpty())
				.distinct()
				.toArray(String[]::new);
	}

}
