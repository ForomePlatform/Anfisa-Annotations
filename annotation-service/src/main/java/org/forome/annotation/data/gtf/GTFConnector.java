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

import org.forome.annotation.config.connector.GTFConfigConnector;
import org.forome.annotation.data.DatabaseConnector;
import org.forome.annotation.data.gtf.struct.GTFRegion;
import org.forome.annotation.data.gtf.struct.GTFResult;
import org.forome.annotation.data.gtf.struct.GTFResultLookup;
import org.forome.annotation.data.gtf.struct.GTFTranscriptRow;
import org.forome.annotation.service.database.DatabaseConnectService;
import org.forome.annotation.utils.DefaultThreadPoolExecutor;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

public class GTFConnector implements AutoCloseable {

    private static final int MAX_THREAD_COUNT = Runtime.getRuntime().availableProcessors();

    private final DatabaseConnector databaseConnector;
    private final GTFDataConnector gtfDataConnector;

    private final ExecutorService threadPoolGTFExecutor;

    public GTFConnector(
            DatabaseConnectService databaseConnectService,
            GTFConfigConnector gtfConfigConnector,
            Thread.UncaughtExceptionHandler uncaughtExceptionHandler
    ) throws Exception {
        this.databaseConnector = new DatabaseConnector(databaseConnectService, gtfConfigConnector);
        this.gtfDataConnector = new GTFDataConnector(databaseConnector);
        threadPoolGTFExecutor = new DefaultThreadPoolExecutor(
                MAX_THREAD_COUNT,
                MAX_THREAD_COUNT,
                0L,
                TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                "GnomadExecutorQueue",
                uncaughtExceptionHandler
        );
    }

    public CompletableFuture<GTFResult> request(String chromosome, long position) {
        CompletableFuture<GTFResult> future = new CompletableFuture();
        threadPoolGTFExecutor.submit(() -> {
            try {
                GTFResult result = gtfDataConnector.getGene(chromosome, position);
                future.complete(result);
            } catch (Throwable e) {
                future.completeExceptionally(e);
            }
        });
        return future;
    }

    public CompletableFuture<GTFRegion> getRegion(String transcript, long position) {
        CompletableFuture<GTFRegion> future = new CompletableFuture();
        threadPoolGTFExecutor.submit(() -> {
            try {
                Object[] result = lookup(position, transcript);
                future.complete((GTFRegion) result[1]);
            } catch (Throwable e) {
                future.completeExceptionally(e);
            }
        });
        return future;
    }

    public CompletableFuture<List<GTFResultLookup>> getRegionByChromosomeAndPositions(String chromosome, long[] positions) {
        CompletableFuture<List<GTFResultLookup>> future = new CompletableFuture();
        threadPoolGTFExecutor.submit(() -> {
            try {
                List<GTFResultLookup> result =lookupByChromosomeAndPositions(chromosome, positions);
                future.complete(result);
            } catch (Throwable e) {
                future.completeExceptionally(e);
            }
        });
        return future;
    }

    public List<GTFTranscriptRow> getTranscriptRows(String transcript) {
        return gtfDataConnector.getTranscriptRows(transcript);
    }

    public Object[] lookup(long pos, String transcript) {
        List<GTFTranscriptRow> rows = gtfDataConnector.getTranscriptRows(transcript);
        if (rows.isEmpty()) return null;

        return lookup(pos, rows);
    }

    public List<GTFResultLookup> lookupByChromosomeAndPositions(String chromosome, long[] positions) {
        List<GTFResultLookup> result = new ArrayList<>();

        List<String> transcripts = gtfDataConnector.getTranscriptsByChromosomeAndPositions(chromosome, positions);
        for (String transcript: transcripts) {
            for (long position: positions) {
                List<GTFTranscriptRow> rows = gtfDataConnector.getTranscriptRows(transcript);
                if (rows.isEmpty()) continue;

                Object[] iResult = lookup(position, transcript);
                GTFRegion region = (GTFRegion)iResult[1];
                result.add(new GTFResultLookup(transcript, rows.get(0).gene, position, region.region, region.indexRegion));
            }
        }

        return result;
    }

    public Object[] lookup(long pos, List<GTFTranscriptRow> rows) {
        long inf = rows.get(0).start;
        if (pos < inf) {
            return new Object[]{(inf - pos), GTFRegion.UPSTREAM};
        }

        long sup = rows.get(rows.size() - 1).end;
        if (pos > sup) {
            return new Object[]{(pos - sup), GTFRegion.DOWNSTREAM};
        }

        List<Integer> a = new ArrayList<>();
        for (GTFTranscriptRow row : rows) {
            a.add(row.start);
            a.add(row.end);
        }

        //Аналог: i = bisect.bisect(a, pos)
        int lo = 0;
        int hi = a.size();
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (pos < a.get(mid)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        int i = lo;

        long d;
        if (pos == inf || pos == sup) {
            d = 0;
        } else {
            d = Math.min(pos - a.get(i - 1), a.get(i) - pos);
        }

        long index;
        String region;
        if ((i % 2) == 1) {
            index = (i + 1) / 2;
            region = "exon";
        } else {
            index = i / 2;
            region = "intron";
        }

        return new Object[]{d, new GTFRegion(region, (int)index), rows.size()};
    }

    @Override
    public void close() {
        databaseConnector.close();
    }
}
