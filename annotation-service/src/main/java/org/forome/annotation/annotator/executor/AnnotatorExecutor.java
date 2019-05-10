package org.forome.annotation.annotator.executor;

import org.forome.annotation.connector.anfisa.AnfisaConnector;
import org.forome.annotation.struct.Sample;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;

public class AnnotatorExecutor implements AutoCloseable {

    private final ThreadExecutor[] threadExecutors;

    private int activeExecutor;

    public AnnotatorExecutor(
            AnfisaConnector anfisaConnector,
            String caseSequence, Map<String, Sample> samples,
            Path pathVepVcf, Path pathVepJson,
            int start, int thread
    ) {
        if (thread < 1) throw new IllegalArgumentException();

        threadExecutors = new ThreadExecutor[thread];
        for (int i = 0; i < thread; i++) {
            threadExecutors[i] = new ThreadExecutor(
                    anfisaConnector,
                    caseSequence, samples,
                    pathVepVcf, pathVepJson,
                    start + i, thread
            );
        }

        activeExecutor = 0;
    }

    public synchronized Result next() {
        ThreadExecutor threadExecutor = threadExecutors[activeExecutor];
        activeExecutor++;
        if (activeExecutor > threadExecutors.length - 1) {
            activeExecutor = 0;
        }

        return threadExecutor.next();
    }

    @Override
    public void close() throws IOException {
        for (ThreadExecutor threadExecutor: threadExecutors) {
            threadExecutor.close();
        }
    }
}
