package org.forome.annotation.annotator;

import io.reactivex.Observable;
import org.forome.annotation.annotator.executor.AnnotatorExecutor;
import org.forome.annotation.annotator.executor.Result;
import org.forome.annotation.annotator.struct.AnnotatorResult;
import org.forome.annotation.annotator.utils.CaseUtils;
import org.forome.annotation.connector.anfisa.AnfisaConnector;
import org.forome.annotation.connector.anfisa.struct.AnfisaResult;
import org.forome.annotation.service.ensemblvep.EnsemblVepService;
import org.forome.annotation.struct.CasePlatform;
import org.forome.annotation.struct.sample.Samples;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;

public class Annotator {

    private final static Logger log = LoggerFactory.getLogger(Annotator.class);

    private static final int MAX_THREAD_COUNT = 1;// Runtime.getRuntime().availableProcessors() * 4;

    private final EnsemblVepService ensemblVepService;
    private final AnfisaConnector anfisaConnector;

    public Annotator(
            EnsemblVepService ensemblVepService,
            AnfisaConnector anfisaConnector) {
        this.ensemblVepService = ensemblVepService;
        this.anfisaConnector = anfisaConnector;
    }

    public AnnotatorResult exec(
            String caseName,
            CasePlatform casePlatform,
            Path pathFam,
            Path pathFamSampleName,
            Path pathVepVcf,
            Path pathVepJson,
            Path cnvFile,
            int startPosition
    ) throws IOException {
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
             InputStream isFamSampleName = (pathFamSampleName != null) ? Files.newInputStream(pathFamSampleName) : null
        ) {
            return exec(
                    caseName,
                    casePlatform,
                    isFam,
                    isFamSampleName,
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
            Path pathVepVcf,
            Path pathVepJson,
            Path cnvFile,
            int startPosition
    ) throws IOException {

        Samples samples = CaseUtils.parseFamFile(isFam, isFamSampleName);

        String caseId = String.format("%s_%s", caseName, casePlatform.name().toLowerCase());

        return annotateJson(
                caseId, samples,
                pathVepVcf, pathVepJson,
                cnvFile,
                startPosition
        );
    }

    public AnnotatorResult annotateJson(
            String caseSequence, Samples samples,
            Path pathVepVcf, Path pathVepJson,
            Path cnvFile,
            int startPosition
    ) {
        return new AnnotatorResult(
                AnnotatorResult.Metadata.build(caseSequence, pathVepVcf, samples, anfisaConnector),
                Observable.create(o -> {
                    new Thread(new Runnable() {
                        @Override
                        public void run() {
                            try (AnnotatorExecutor annotatorExecutor = new AnnotatorExecutor(
                                    ensemblVepService, anfisaConnector,
                                    caseSequence, samples,
                                    pathVepVcf, pathVepJson,
                                    cnvFile,
                                    startPosition, MAX_THREAD_COUNT,
                                    (t, e) -> o.tryOnError(e)
                            )) {
                                boolean run = true;
                                while (run) {
                                    Result result = annotatorExecutor.next();
                                    AnfisaResult anfisaResult;
                                    try {
                                        anfisaResult = result.future.get();
                                        if (anfisaResult != null) {
                                            o.onNext(anfisaResult);
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
