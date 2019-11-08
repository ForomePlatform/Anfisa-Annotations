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

package org.forome.annotation.service.ensemblvep.external;

import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.HttpRequest;
import org.apache.http.HttpResponse;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.concurrent.FutureCallback;
import org.apache.http.impl.nio.client.CloseableHttpAsyncClient;
import org.apache.http.impl.nio.client.HttpAsyncClients;
import org.apache.http.impl.nio.conn.PoolingNHttpClientConnectionManager;
import org.apache.http.impl.nio.reactor.DefaultConnectingIOReactor;
import org.apache.http.message.BasicHeader;
import org.apache.http.util.EntityUtils;
import org.forome.annotation.exception.ExceptionBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Closeable;
import java.io.IOException;
import java.net.URI;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentLinkedDeque;

public class EnsemblVepHttpClient implements Closeable {

	private final static Logger log = LoggerFactory.getLogger(EnsemblVepHttpClient.class);

	private class QueueRequest {

		public final String endpoint;
		public final CompletableFuture<JSONArray> future;

		public QueueRequest(String endpoint, CompletableFuture<JSONArray> future) {
			this.endpoint = endpoint;
			this.future = future;
		}
	}

	private static String HOST = "grch37.rest.ensembl.org";

	private static int MAX_REQUEST_IN_SECOND = 15;
	private static long REQUEST_PAUSE = 1000L / MAX_REQUEST_IN_SECOND;

	private final RequestConfig requestConfig;
	private final PoolingNHttpClientConnectionManager connectionManager;

	private final HttpHost httpHost;

	private final ConcurrentLinkedDeque<QueueRequest> qequeRequests;

	private boolean active = true;

	private Thread.UncaughtExceptionHandler uncaughtExceptionHandler;

	protected EnsemblVepHttpClient(
			Thread.UncaughtExceptionHandler uncaughtExceptionHandler
	) throws IOException {
		requestConfig = RequestConfig.custom()
				.setConnectTimeout(2000)//Таймаут на подключение
				.setSocketTimeout(1 * 60 * 1000)//Таймаут между пакетами
				.setConnectionRequestTimeout(1 * 60 * 1000)//Таймаут на ответ
				.build();

		connectionManager = new PoolingNHttpClientConnectionManager(new DefaultConnectingIOReactor());
		connectionManager.setMaxTotal(100);
		connectionManager.setDefaultMaxPerRoute(100);

		httpHost = new HttpHost(HOST);

		qequeRequests = new ConcurrentLinkedDeque<>();

		this.uncaughtExceptionHandler = uncaughtExceptionHandler;

		Thread thread = new Thread(new Runnable() {
			@Override
			public void run() {
				while (active) {
					execute();
					try {
						Thread.sleep(REQUEST_PAUSE);
					} catch (InterruptedException e) {
					}
				}
			}
		});
		thread.setDaemon(true);
		thread.start();
	}

	protected CompletableFuture<JSONArray> request(String endpoint) {
		CompletableFuture<JSONArray> future = new CompletableFuture<>();
		qequeRequests.add(new QueueRequest(endpoint, future));
		return future;
	}

	private void execute() {
		QueueRequest queueRequest = qequeRequests.poll();
		if (queueRequest == null) return;

		String endpoint = queueRequest.endpoint;
		CompletableFuture<JSONArray> future = queueRequest.future;

		try {

			HttpRequest httpRequest = new HttpGet(new URI("http://" + HOST + endpoint));
			httpRequest.addHeader(new BasicHeader("Content-Type", "application/json"));

			CloseableHttpAsyncClient httpclient = HttpAsyncClients.custom()
					.setDefaultRequestConfig(requestConfig)
//					.setConnectionManager(connectionManager)
//					.setConnectionManagerShared(true)
					.build();
			httpclient.start();

			httpclient.execute(httpHost, httpRequest, new FutureCallback<HttpResponse>() {
				@Override
				public void completed(HttpResponse response) {
					try {
						HttpEntity entity = response.getEntity();
						String entityBody = EntityUtils.toString(entity);

						Object rawResponse;
						try {
							rawResponse = new JSONParser(JSONParser.DEFAULT_PERMISSIVE_MODE).parse(entityBody);
						} catch (Exception e) {
							throw ExceptionBuilder.buildExternalServiceException(e, "Error parse response, endpoint: " + endpoint + " response: '" + entityBody + "'");
						}

						if (rawResponse instanceof JSONArray) {
							future.complete((JSONArray) rawResponse);
						} else if (rawResponse instanceof JSONObject && ((JSONObject) rawResponse).containsKey("error")) {
							throw ExceptionBuilder.buildExternalServiceException(new RuntimeException(), "Error parse response, endpoint: " + endpoint + " response: '" + entityBody + "'");
						} else if (response.getStatusLine().getStatusCode() == 429) {
							//Повторно кладем в очередь на расчет
							log.warn("External service response: Too many request - to repeat: {}", endpoint);
							qequeRequests.addFirst(queueRequest);
						} else {
							uncaughtExceptionHandler.uncaughtException(
									Thread.currentThread(),
									new IOException("Unknown response, endpoint: " + endpoint + " response: '" + entityBody + "'")
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
					log.warn("External service: exception {} - to repeat: {}", ex.getMessage(), endpoint);
					qequeRequests.addFirst(queueRequest);
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
		} catch (Throwable ex) {
			log.error("Exception close connect", ex);
			future.completeExceptionally(ex);
		}
	}

	@Override
	public void close() {
		active = false;
	}
}
