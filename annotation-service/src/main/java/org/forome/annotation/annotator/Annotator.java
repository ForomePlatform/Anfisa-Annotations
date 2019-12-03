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

package org.forome.annotation.annotator;

import io.reactivex.Observable;
import net.minidev.json.parser.ParseException;
import org.forome.annotation.annotator.executor.AnnotatorExecutor;
import org.forome.annotation.annotator.executor.Result;
import org.forome.annotation.annotator.struct.AnnotatorResult;
import org.forome.annotation.annotator.utils.CaseUtils;
import org.forome.annotation.processing.Processing;
import org.forome.annotation.processing.struct.ProcessingResult;
import org.forome.annotation.service.ensemblvep.EnsemblVepService;
import org.forome.annotation.struct.CasePlatform;
import org.forome.annotation.struct.mcase.MCase;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class Annotator {

	private final static Logger log = LoggerFactory.getLogger(Annotator.class);

	private static final int MAX_THREAD_COUNT = Runtime.getRuntime().availableProcessors() * 4;

	private final EnsemblVepService ensemblVepService;
	private final Processing processing;

	public Annotator(
			EnsemblVepService ensemblVepService,
			Processing processing) {
		this.ensemblVepService = ensemblVepService;
		this.processing = processing;
	}

	public AnnotatorResult exec(
			String caseName,
			CasePlatform casePlatform,
			Path pathFam,
			Path pathFamSampleName,
			Path pathCohorts,
			Path pathVepVcf,
			Path pathVepJson,
			Path cnvFile,
			int startPosition
	) throws IOException, ParseException {
		if (!Files.exists(pathFam)) {
			throw new RuntimeException("Fam file is not exists: " + pathFam.toAbsolutePath());
		}
		if (!pathFam.getFileName().toString().endsWith(".fam")) {
			throw new IllegalArgumentException("Bad name fam file: " + pathFam.toAbsolutePath());
		}

		if (!Files.exists(pathVepVcf)) {
			throw new RuntimeException("Vcf file is not exists: " + pathVepVcf.toAbsolutePath());
		}
		if (!pathVepVcf.getFileName().toString().endsWith(".vcf")) {
			throw new IllegalArgumentException("Bad name vcf file (Need *.vcf): " + pathVepVcf.toAbsolutePath());
		}

		if (pathVepJson != null) {
			if (!Files.exists(pathVepJson)) {
				throw new RuntimeException("VepJson file is not exists: " + pathVepJson.toAbsolutePath());
			}
			String vepJsonFileName = pathVepJson.getFileName().toString();
			if (!(vepJsonFileName.endsWith(".json") || vepJsonFileName.endsWith(".json.gz"))) {
				throw new IllegalArgumentException("Bad name pathVepJson file (Need *.json vs *.json.gz): " + pathVepJson.toAbsolutePath());
			}
		}

		try (InputStream isFam = Files.newInputStream(pathFam);
			 InputStream isFamSampleName = (pathFamSampleName != null) ? Files.newInputStream(pathFamSampleName) : null;
			 InputStream isCohorts = (pathCohorts != null) ? Files.newInputStream(pathCohorts) : null
		) {
			return exec(
					caseName,
					casePlatform,
					isFam,
					isFamSampleName,
					isCohorts,
					pathVepVcf,
					pathVepJson,
					cnvFile,
					startPosition
			);
		}
	}

	public AnnotatorResult exec(
			String caseName,
			CasePlatform casePlatform,
			InputStream isFam,
			InputStream isFamSampleName,
			InputStream isCohorts,
			Path pathVepVcf,
			Path pathVepJson,
			Path cnvFile,
			int startPosition
	) throws IOException, ParseException {

		MCase samples = CaseUtils.parseFamFile(isFam, isFamSampleName, isCohorts);

		String caseId = String.format("%s_%s", caseName, casePlatform.name().toLowerCase());

		return annotateJson(
				caseId, samples,
				pathVepVcf, pathVepJson,
				cnvFile,
				startPosition
		);
	}

	public AnnotatorResult annotateJson(
			String caseSequence, MCase samples,
			Path pathVepVcf, Path pathVepJson,
			Path cnvFile,
			int startPosition
	) {
		return new AnnotatorResult(
				AnnotatorResult.Metadata.build(caseSequence, pathVepVcf, samples, processing.getAnfisaConnector()),
				Observable.create(o -> {
					new Thread(new Runnable() {
						@Override
						public void run() {
							try (AnnotatorExecutor annotatorExecutor = new AnnotatorExecutor(
									ensemblVepService, processing,
									caseSequence, samples,
									pathVepVcf, pathVepJson,
									cnvFile,
									startPosition, MAX_THREAD_COUNT,
									(t, e) -> o.tryOnError(e)
							)) {
								boolean run = true;
								while (run) {
									Result result = annotatorExecutor.next();
									List<ProcessingResult> processingResults;
									try {
										processingResults = result.future.get();
										if (processingResults != null) {
											for (ProcessingResult processingResult : processingResults) {
												o.onNext(processingResult);
											}
										} else {
											run = false;
										}
									} catch (Throwable e) {
										log.error("throwable", e);
									}
								}
								o.onComplete();
							} catch (Throwable e) {
								o.tryOnError(e);
							}
						}
					}).start();
				})
		);
	}

}
