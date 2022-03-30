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

package org.forome.annotation.service.source.external.fasta;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.HttpResponse;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.concurrent.FutureCallback;
import org.apache.http.impl.nio.client.CloseableHttpAsyncClient;
import org.apache.http.impl.nio.client.HttpAsyncClients;
import org.apache.http.impl.nio.conn.PoolingNHttpClientConnectionManager;
import org.apache.http.impl.nio.reactor.DefaultConnectingIOReactor;
import org.apache.http.nio.reactor.IOReactorException;
import org.apache.http.util.EntityUtils;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.utils.Statistics;
import org.forome.core.struct.Assembly;
import org.forome.core.struct.Interval;
import org.forome.core.struct.sequence.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.net.URI;
import java.net.URL;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

/**
 * curl "localhost:8290/get?array=fasta&type=hg19&loc=2:73675217-73675248"
 * curl "localhost:8290/get?array=fasta&type=hg38&loc=2:73448087-73448121"
 */
public class FastaHttpRequest {

	private final static Logger log = LoggerFactory.getLogger(FastaHttpRequest.class);

	private final URL url;

	private final RequestConfig requestConfig;
	private final PoolingNHttpClientConnectionManager connectionManager;
	private final HttpHost httpHost;

	private final Cache cache;

	private final Statistics statistics;

	public FastaHttpRequest(
			URL url
	) throws IOReactorException {
		this.url = url;

		requestConfig = RequestConfig.custom()
				.setConnectTimeout(5000)//Таймаут на подключение
				.setSocketTimeout(10 * 60 * 1000)//Таймаут между пакетами
				.setConnectionRequestTimeout(10 * 60 * 1000)//Таймаут на ответ
				.build();

		connectionManager = new PoolingNHttpClientConnectionManager(new DefaultConnectingIOReactor());
		connectionManager.setMaxTotal(100);
		connectionManager.setDefaultMaxPerRoute(100);

		httpHost = new HttpHost(url.getHost(), url.getPort(), "http");

		this.cache = CacheBuilder.newBuilder()
				.maximumSize(1000)
				.build();

		this.statistics = new Statistics();
	}

	public Sequence getSequence(Assembly assembly, Interval interval) {
		String key = String.format(
				"%s:%s:%s:%s", assembly.name(), interval.chromosome.getChar(),
				interval.start, interval.end
		);

		long t1 = System.currentTimeMillis();
		String value;
		try {
			value = ((Optional<String>) cache.get(key, () -> {
				JSONObject response = request(
						String.format("http://%s:%s/get?array=fasta&type=%s&loc=%s:%s-%s",
								url.getHost(), url.getPort(),
								(assembly == Assembly.GRCh37) ? "hg19" : "hg38",
								interval.chromosome.getChar(), interval.start, interval.end
						)
				);
				return Optional.ofNullable(response.getAsString("fasta"));
			})).orElse(null);
		} catch (ExecutionException e) {
			throw new RuntimeException(e);
		} finally {
			statistics.addTime(System.currentTimeMillis() - t1);
		}

		if (value == null) {
			return null;
		} else {
			return Sequence.build(interval, value);
		}
	}

	private JSONObject request(String url) {
		CompletableFuture<JSONObject> future = new CompletableFuture<>();
		try {
			CloseableHttpAsyncClient httpclient = HttpAsyncClients.custom()
					.setDefaultRequestConfig(requestConfig)
					.build();
			httpclient.start();

			HttpPost httpPostRequest = new HttpPost(new URI(url));

			httpclient.execute(httpHost, httpPostRequest, new FutureCallback<HttpResponse>() {
				@Override
				public void completed(HttpResponse response) {
					try {
						HttpEntity entity = response.getEntity();
						String entityBody = EntityUtils.toString(entity);

						Object rawResponse;
						try {
							rawResponse = new JSONParser(JSONParser.DEFAULT_PERMISSIVE_MODE).parse(entityBody);
						} catch (Exception e) {
							throw ExceptionBuilder.buildExternalServiceException(new RuntimeException("Exception parse response external service, response: " + entityBody));
						}
						if (rawResponse instanceof JSONObject) {
							future.complete((JSONObject) rawResponse);
						} else {
							throw ExceptionBuilder.buildExternalServiceException(
									new RuntimeException("Exception external service(AStorage), response: " + entityBody),
									"AStorage", "Response: " + entityBody
							);
						}
					} catch (Throwable ex) {
						future.completeExceptionally(ex);
					}

					try {
						httpclient.close();
					} catch (IOException ignore) {
						log.error("Exception close connect");
					}
				}

				@Override
				public void failed(Exception ex) {
					future.completeExceptionally(ex);
					try {
						httpclient.close();
					} catch (IOException ignore) {
						log.error("Exception close connect");
					}
				}

				@Override
				public void cancelled() {
					future.cancel(true);
					try {
						httpclient.close();
					} catch (IOException ignore) {
						log.error("Exception close connect");
					}
				}
			});
		} catch (Throwable t) {
			log.error("Exception close connect", t);
			future.completeExceptionally(t);
		}

		try {
			return future.get();
		} catch (InterruptedException e) {
			throw new RuntimeException();
		} catch (ExecutionException e) {
			String msg = "Exception connecting to " + url + ": " + e.getMessage ();
			log.error (msg);
			throw new RuntimeException(msg, e.getCause());
		}
	}

	public Statistics.Stat getStatistics() {
		return statistics.getStat();
	}
}
