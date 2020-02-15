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

package org.forome.annotation.makedatabase.makesourcedata.conservation;

import org.apache.commons.lang3.RandomUtils;
import org.forome.annotation.service.database.struct.batch.BatchRecord;
import org.forome.annotation.service.database.struct.batch.BatchRecordConservation;
import org.forome.annotation.struct.Chromosome;
import org.forome.annotation.struct.Interval;
import org.forome.annotation.struct.Position;
import org.junit.Assert;
import org.junit.Test;

import java.util.Objects;

public class MakeConservationBuildTest {

	@Test
	public void test() {
		for (int k = 0; k < 100; k += 23) {

			Interval interval = Interval.of(
					Chromosome.CHR_1,
					k * BatchRecord.DEFAULT_SIZE, (k + 1) * BatchRecord.DEFAULT_SIZE - 1
			);

			for (int t = 0; t < 10000; t++) {

				//generate
				MakeConservationBuild.Data[] values = new MakeConservationBuild.Data[BatchRecord.DEFAULT_SIZE];
				for (int i = 0; i < values.length; i++) {
					MakeConservationBuild.Data data = new MakeConservationBuild.Data();
					data.priPhCons = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.mamPhCons = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.verPhCons = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.priPhyloP = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.mamPhyloP = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.verPhyloP = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.gerpRS = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.gerpRSpval = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.gerpN = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					data.gerpS = (RandomUtils.nextBoolean())?null:RandomUtils.nextFloat(0.0f, 62.0f) - 31.0f;
					values[i] = data;
				}
				MakeConservationBuild makeConservationBuild = new MakeConservationBuild(interval, values);
				byte[] bytes = makeConservationBuild.build();

				//restore
				BatchRecordConservation batchRecordConservation = new BatchRecordConservation(interval, bytes, 0);

				//assert
				for (int p = interval.start; p < interval.end; p++) {
					Position position = new Position(interval.chromosome, p);

					assertFloat(values[p - interval.start].priPhCons, batchRecordConservation.getPriPhCons(position));
					assertFloat(values[p - interval.start].mamPhCons, batchRecordConservation.getMamPhCons(position));
					assertFloat(values[p - interval.start].verPhCons, batchRecordConservation.getVerPhCons(position));
					assertFloat(values[p - interval.start].priPhyloP, batchRecordConservation.getPriPhyloP(position));
					assertFloat(values[p - interval.start].mamPhyloP, batchRecordConservation.getMamPhyloP(position));
					assertFloat(values[p - interval.start].verPhyloP, batchRecordConservation.getVerPhyloP(position));
					assertFloat(values[p - interval.start].gerpRS, batchRecordConservation.getGerpRS(position));
					assertFloat(values[p - interval.start].gerpRSpval, batchRecordConservation.getGerpRSpval(position));
					assertFloat(values[p - interval.start].gerpN, batchRecordConservation.getGerpN(position));
					assertFloat(values[p - interval.start].gerpS, batchRecordConservation.getGerpS(position));
				}
			}
		}
	}

	private void assertFloat(Float expected, Float actual) {
		if (Objects.equals(expected, actual)) return;
		Assert.assertEquals(expected, actual, 0.001d);
	}
}