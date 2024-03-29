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

package org.forome.annotation.data.gtf;

import org.forome.annotation.data.gtf.mysql.struct.GTFRegion;
import org.forome.annotation.data.gtf.mysql.struct.GTFResultLookup;
import org.forome.core.struct.Assembly;
import org.forome.core.struct.Chromosome;
import org.forome.core.struct.Position;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class GTFSimpleTest extends GTFBaseTest {

	@Test
	public void testByChromosomeAndPositions() throws Exception {
		String chromosome = "5";
		String transcript = "ENST00000282356";
		int position = 110694251;

		GTFRegion expectedGtfRegion = gtfConnector.getRegion(
				null,
				Assembly.GRCh37,
				new Position(Chromosome.of(chromosome), position),
				transcript
		).get();


		List<GTFResultLookup> lookups = gtfConnector
				.getRegionByChromosomeAndPositions(null, chromosome, new long[] {position}).get();
		GTFResultLookup actualGTFResultLookup = lookups.stream()
				.filter(gtfResultLookup -> transcript.equals(gtfResultLookup.transcript)).findFirst().get();

		Assert.assertEquals(expectedGtfRegion.region, actualGTFResultLookup.region);
		Assert.assertEquals(expectedGtfRegion.indexRegion, actualGTFResultLookup.index);
	}
}
