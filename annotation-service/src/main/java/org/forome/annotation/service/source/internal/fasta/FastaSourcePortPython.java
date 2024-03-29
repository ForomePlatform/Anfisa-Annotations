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

package org.forome.annotation.service.source.internal.fasta;

import org.forome.astorage.pastorage.PAStorage;
import org.forome.astorage.pastorage.record.RecordFasta;
import org.forome.astorage.pastorage.schema.SchemaFasta;
import org.forome.core.struct.Assembly;
import org.forome.core.struct.Interval;
import org.forome.core.struct.Position;
import org.forome.core.struct.nucleotide.Nucleotide;
import org.forome.core.struct.sequence.Sequence;

public class FastaSourcePortPython {

	private final PAStorage paStorage;

	public FastaSourcePortPython(PAStorage paStorage) {
		this.paStorage = paStorage;
	}

	public Sequence getSequence(Assembly assembly, Interval interval) {
		SchemaFasta schemaFasta = (SchemaFasta) paStorage.getSchema(SchemaFasta.SCHEMA_FASTA_NAME);

		Nucleotide[] nucleotides = new Nucleotide[interval.end - interval.start + 1];
		for (int i = 0; i < nucleotides.length; i++) {
			Position position = new Position(interval.chromosome, interval.start + i);
			RecordFasta recordFasta = schemaFasta.getRecord(assembly, position);

			Nucleotide nucleotide;
			if (recordFasta == null) {
				nucleotide = Nucleotide.NONE;
			} else {
				nucleotide = recordFasta.nucleotide;
				if (nucleotide == null) {
					nucleotide = Nucleotide.NONE;
				}
			}
			nucleotides[i] = nucleotide;
		}

		Sequence sequence = new Sequence(interval, nucleotides);

		return sequence;
	}
}
