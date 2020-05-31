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

package org.forome.annotation.data.gnomad;

import com.google.common.collect.ImmutableList;
import org.forome.annotation.data.gnomad.datasource.GnomadDataSource;
import org.forome.annotation.data.gnomad.struct.DataResponse;
import org.forome.annotation.data.gnomad.struct.GnamadGroup;
import org.forome.annotation.data.gnomad.struct.GnomadResult;
import org.forome.annotation.matcher.SequenceMatcher;
import org.forome.annotation.struct.Assembly;
import org.forome.annotation.struct.Chromosome;
import org.forome.annotation.struct.SourceMetadata;
import org.forome.annotation.utils.DefaultThreadPoolExecutor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

public class GnomadConnectorImpl implements GnomadConnector {

	private final static Logger log = LoggerFactory.getLogger(GnomadConnectorImpl.class);

	private static final int MAX_THREAD_COUNT = Runtime.getRuntime().availableProcessors();

	private final GnomadDataSource gnomadDataSource;

	private final ExecutorService threadPoolGnomadExecutor;

	public GnomadConnectorImpl(
			GnomadDataSource gnomadDataSource,
			Thread.UncaughtExceptionHandler uncaughtExceptionHandler
	) throws Exception {
		this.gnomadDataSource = gnomadDataSource;
//		databaseConnector = new DatabaseConnector(databaseConnectService, gnomadConfigConnector);
//		gnomadDataConnector = new GnomadDataConnector(databaseConnector);
		threadPoolGnomadExecutor = new DefaultThreadPoolExecutor(
				MAX_THREAD_COUNT,
				MAX_THREAD_COUNT,
				0L,
				TimeUnit.MILLISECONDS,
				new LinkedBlockingQueue<>(),
				"GnomadExecutorQueue",
				uncaughtExceptionHandler
		);
	}

	public CompletableFuture<GnomadResult> request(Assembly assembly, Chromosome chromosome, int position, String reference, String alternative) {
		CompletableFuture<GnomadResult> future = new CompletableFuture();
		threadPoolGnomadExecutor.submit(() -> {
			try {
				GnomadResult result = syncRequest(assembly, chromosome, position, reference, alternative);
				future.complete(result);
			} catch (Throwable e) {
				future.completeExceptionally(e);
			}
		});
		return future;
	}

	public List<SourceMetadata> getSourceMetadata() {
		return gnomadDataSource.getSourceMetadata();
	}

	private GnomadResult syncRequest(Assembly assembly, Chromosome chromosome, int position, String reference, String alternative) throws Exception {
		List<DataResponse> exomes = gnomadDataSource.getData(
				assembly, chromosome, position, reference, alternative, "e"
		);
		List<DataResponse> genomes = gnomadDataSource.getData(
				assembly, chromosome, position, reference, alternative, "g"
		);

		List<DataResponse> overall = new ImmutableList.Builder().addAll(exomes).addAll(genomes).build();
		if (overall.isEmpty()) {
			return null;
		}

		GnomadResult.Sum sumExomes = null;
		if (!exomes.isEmpty()) {
			long an = countAN(exomes, null);
			long ac = countAC(exomes, null);
			double af = countAF(an, ac);
			long hom = countHom(exomes);
			Long hem = countHem(exomes);
			sumExomes = new GnomadResult.Sum(an, ac, af, hom, hem);
		}

		GnomadResult.Sum sumGenomes = null;
		if (!genomes.isEmpty()) {
			long an = countAN(genomes, null);
			long ac = countAC(genomes, null);
			double af = countAF(an, ac);
			long hom = countHom(genomes);
			Long hem = countHem(genomes);
			sumGenomes = new GnomadResult.Sum(an, ac, af, hom, hem);
		}

		GnomadResult.Sum sumOverall = null;
		if (!overall.isEmpty()) {
			long an = countAN(overall, null);
			long ac = countAC(overall, null);
			double af = countAF(an, ac);
			long hom = countHom(overall);
			Long hem = countHem(overall);
			sumOverall = new GnomadResult.Sum(an, ac, af, hom, hem);
		}

		GnomadResult.Popmax popmax = countPopmaxFromRows(overall, GnamadGroup.Type.GENERAL);
		GnomadResult.Popmax widePopmax = countPopmaxFromRows(overall, null);

		Set<GnomadResult.Url> urls = new HashSet<>();
		for (DataResponse item : overall) {
			String chrom = item.getValue("CHROM");
			long pos = ((Number) item.getValue("POS")).longValue();
			String ref = item.getValue("REF");
			String alt = item.getValue("ALT");

			SequenceMatcher matcher = new SequenceMatcher(ref, alt);
			List<SequenceMatcher.Tuple3<Integer, Integer, Integer>> matches = matcher.getMatchingBlocks();
			List<String> a = new ArrayList<>();
			List<String> b = new ArrayList<>();
			int apos = 0;
			int bpos = 0;
			for (SequenceMatcher.Tuple3<Integer, Integer, Integer> match : matches) {
				a.add(ref.substring(apos, match.value0));
				apos = match.value0 + match.value2;
				b.add(alt.substring(bpos, match.value1));
				bpos = match.value1 + match.value2;
			}

			String newRef = String.join("", a);
			String newAlt = String.join("", b);

			if (newRef.isEmpty() || newAlt.isEmpty()) {
				String firstMatch = String.valueOf(ref.charAt(0));
				newRef = firstMatch + newRef;
				newAlt = firstMatch + newAlt;
			}

			urls.add(new GnomadResult.Url(chrom, pos, newRef, newAlt));
		}

		return new GnomadResult(
				sumExomes, sumGenomes, sumOverall,
				popmax, widePopmax,
				urls
		);
	}

	@Override
	public void close() {
		gnomadDataSource.close();
	}

	private static long countAN(List<DataResponse> items, GnamadGroup group) {
		long an = 0;
		String anColumn;
		if (group == null) {
			anColumn = "AN";
		} else {
			anColumn = "AN_" + group.name();
		}
		for (DataResponse item : items) {
			Number value = item.getValue(anColumn);
			if (value == null) continue;
			an += value.longValue();
		}
		return an;
	}

	private static long countAC(List<DataResponse> items, GnamadGroup group) {
		long ac = 0;
		String acColumn;
		if (group == null) {
			acColumn = "AC";
		} else {
			acColumn = "AC_" + group;
		}
		for (DataResponse item : items) {
			Number value = item.getValue(acColumn);
			if (value == null) continue;
			ac += value.longValue();
		}
		return ac;
	}

	private static double countAF(long an, long ac) {
		double af;
		if (an > 0) {
			af = (double) ac / (double) an;
		} else {
			af = 0;
		}
		return af;
	}

	private static GnomadResult.Popmax countPopmaxFromRows(List<DataResponse> overall, GnamadGroup.Type type) {
		GnamadGroup group = null;
		Double popmaxAF = null;
		long popmaxAN = 0;

		GnamadGroup[] groups = (type == null) ? GnamadGroup.values() : GnamadGroup.getByType(type);
		for (GnamadGroup iGroup : groups) {
			long an = countAN(overall, iGroup);
			long ac = countAC(overall, iGroup);
			if (an == 0) {
				continue;
			}
			double af = (double) ac / ((double) an);
			if (popmaxAF == null || af > popmaxAF) {
				group = iGroup;
				popmaxAF = af;
				popmaxAN = an;
			}
		}

		if (group == null) {
			return null;
		} else {
			return new GnomadResult.Popmax(
					group, popmaxAF, popmaxAN
			);
		}
	}

	private static long countHom(List<DataResponse> items) {
		long hom = 0;
		String column = "nhomalt";
		for (DataResponse item : items) {
			Number value = item.getValue(column);
			if (value == null) continue;
			hom += value.longValue();
		}
		return hom;
	}

	private static Long countHem(List<DataResponse> items) {
		Long hem = null;
		String column = "hem";
		for (DataResponse item : items) {
			Number value = item.getValue(column);
			if (value == null) continue;
			if (hem == null) {
				hem = value.longValue();
			} else {
				hem += value.longValue();
			}
		}
		return hem;
	}

//    private static Long countHem(String chromosome, List<GnomadDataConnector.Result> items) {
//        if ("X".equals(chromosome.toUpperCase())) {
//            return countAC(items, "Male");
//        }
//        return null;
//    }
}
