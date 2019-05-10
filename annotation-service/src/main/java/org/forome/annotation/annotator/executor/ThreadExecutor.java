package org.forome.annotation.annotator.executor;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.minidev.json.JSONObject;
import org.forome.annotation.annotator.input.VepJsonIterator;
import org.forome.annotation.connector.anfisa.AnfisaConnector;
import org.forome.annotation.connector.anfisa.struct.AnfisaResult;
import org.forome.annotation.controller.utils.RequestParser;
import org.forome.annotation.struct.Sample;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.Deque;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentLinkedDeque;

public class ThreadExecutor implements AutoCloseable {

    private final AnfisaConnector anfisaConnector;

    private final String caseSequence;
    private final Map<String, Sample> samples;

    private final Path pathVepVcf;
    private final Path pathVepJson;

    private final int start;
    private final int step;

    private final VCFFileReader vcfFileReader;
    private final CloseableIterator<VariantContext> vcfFileReaderIterator;

    private final VepJsonIterator vepJsonIterator;

    private Result nextResult;
    private final Deque<Result> waitExecuteVariants;//Варианты ожидающие выполнения

    private int nextPosition;
    private boolean isCompleted = false;

    public ThreadExecutor(
            AnfisaConnector anfisaConnector,
            String caseSequence, Map<String, Sample> samples,
            Path pathVepVcf, Path pathVepJson,
            int start, int step,
            Thread.UncaughtExceptionHandler uncaughtExceptionHandler
    ) {
        this.anfisaConnector = anfisaConnector;

        this.caseSequence = caseSequence;
        this.samples = samples;

        this.pathVepVcf = pathVepVcf;
        this.pathVepJson = pathVepJson;

        this.start = start;
        this.step = step;

        this.vcfFileReader = new VCFFileReader(pathVepVcf, false);
        this.vcfFileReaderIterator = vcfFileReader.iterator();

        if (pathVepJson != null) {
            vepJsonIterator = new VepJsonIterator(pathVepJson);
        } else {
            vepJsonIterator = null;
        }

        this.nextResult = new Result(nextPosition, new CompletableFuture<>());
        this.waitExecuteVariants = new ConcurrentLinkedDeque<>();
        this.waitExecuteVariants.add(nextResult);

        nextPosition = start + step;

        //Исполнитель
        Thread executor = new Thread(() -> {
            //Прокручиваем до начала итерации
            Source source;
            if (start > 0) {
                nextSource(start);
            }
            source = nextSource(1);

            while (!isCompleted) {

                //TODO Переписать на засыпание потока
                while (waitExecuteVariants.isEmpty()) {
                    if (isCompleted) return;
                    try {
                        Thread.sleep(10L);
                    } catch (InterruptedException e) {
                    }
                }
                Result result = waitExecuteVariants.poll();
                if (isCompleted) {
                    result.future.complete(null);
                    continue;
                }

                VariantContext variantContext = source.variantContext;
                JSONObject vepJson = source.getVepJson();

                if (vepJsonIterator != null) {
                    String chromosome = RequestParser.toChromosome(vepJson.getAsString("seq_region_name"));
                    long iStart = vepJson.getAsNumber("start").longValue();
                    long iEnd = vepJson.getAsNumber("end").longValue();

                    AnfisaResult anfisaResult = anfisaConnector.build(caseSequence, chromosome, iStart, iEnd, vepJson, variantContext, samples);
                    result.future.complete(anfisaResult);
                } else {
                    String chromosome = RequestParser.toChromosome(variantContext.getContig());
                    long iStart = variantContext.getStart();
                    long iEnd = variantContext.getEnd();

                    //variantContext.getAltAlleleWithHighestAlleleCount();
                    Allele allele = variantContext.getAlternateAlleles().stream()
                            .filter(iAllele -> !iAllele.getDisplayString().equals("*"))
                            .max(Comparator.comparing(variantContext::getCalledChrCount))
                            .orElse(null);
                    String alternative = allele.getDisplayString();

                    anfisaConnector.request(chromosome, iStart, iEnd, alternative)
                            .thenApply(anfisaResults -> {
                                result.future.complete(anfisaResults.get(0));
                                return null;
                            })
                            .exceptionally(throwable -> {
                                result.future.completeExceptionally(throwable);
                                return null;
                            });
                }

                source = nextSource(step);
            }
        });
        executor.setUncaughtExceptionHandler(uncaughtExceptionHandler);
        executor.start();
    }

    private Source nextSource(int step) {
        if (step < 1) throw new IllegalArgumentException();
        VariantContext variantContext = null;
        String strVepJson = null;
        for (int i = 0; i < step; i++) {
            variantContext = vcfFileReaderIterator.next();
            strVepJson = (vepJsonIterator != null) ? vepJsonIterator.next() : null;
        }
        return new Source(variantContext, strVepJson);
    }

    public Result next() {
        Result value;
        synchronized (waitExecuteVariants) {
            value = nextResult;

            nextPosition += step;
            if (isCompleted) {
                return new Result(nextPosition, CompletableFuture.completedFuture(null));
            } else {
                nextResult = new Result(nextPosition, new CompletableFuture<>());
                this.waitExecuteVariants.add(nextResult);
            }
        }
        return value;
    }


    @Override
    public void close() throws IOException {
        this.vcfFileReaderIterator.close();
        this.vcfFileReader.close();

        if (vepJsonIterator != null) {
            vepJsonIterator.close();
        }
    }
}
