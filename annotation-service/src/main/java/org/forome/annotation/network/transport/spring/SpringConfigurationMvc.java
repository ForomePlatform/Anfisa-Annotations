package org.forome.annotation.network.transport.spring;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;
import org.springframework.http.HttpStatus;
import org.springframework.util.AntPathMatcher;
import org.springframework.web.bind.annotation.ResponseStatus;
import org.springframework.web.context.request.NativeWebRequest;
import org.springframework.web.context.request.async.DeferredResult;
import org.springframework.web.context.request.async.DeferredResultProcessingInterceptorAdapter;
import org.springframework.web.multipart.commons.CommonsMultipartResolver;
import org.springframework.web.servlet.config.annotation.AsyncSupportConfigurer;
import org.springframework.web.servlet.config.annotation.EnableWebMvc;
import org.springframework.web.servlet.config.annotation.PathMatchConfigurer;
import org.springframework.web.servlet.config.annotation.WebMvcConfigurerAdapter;

import java.time.Duration;

/**
 * Created by kris on 13.10.16.
 */
@EnableWebMvc
@Configuration
@ComponentScan({ "org.forome.annotation" })
public class SpringConfigurationMvc extends WebMvcConfigurerAdapter {

	private static Duration requestTimeout;

	@Override
	public void configureAsyncSupport(AsyncSupportConfigurer configurer) {
		configurer.setDefaultTimeout(requestTimeout.toMillis());
		configurer.registerDeferredResultInterceptors(
				new DeferredResultProcessingInterceptorAdapter() {
					@Override
					public <T> boolean handleTimeout(NativeWebRequest req, DeferredResult<T> result) {
						return result.setErrorResult(new AsyncTimeoutException());
					}
				});
	}

	@Override
	public void configurePathMatch(PathMatchConfigurer configurer) {
		AntPathMatcher matcher = new AntPathMatcher();
		matcher.setCaseSensitive(false);
		configurer.setPathMatcher(matcher);
	}

	@Bean
	public CommonsMultipartResolver multipartResolver() {
		CommonsMultipartResolver commonsMultipartResolver = new CommonsMultipartResolver();
		commonsMultipartResolver.setDefaultEncoding("utf-8");
		commonsMultipartResolver.setMaxUploadSize(-1);//Убираем ограничение на размер загружаемых файлов
		return commonsMultipartResolver;
	}

	public static void init(Duration requestTimeout) {
		SpringConfigurationMvc.requestTimeout = requestTimeout;
	}

	@ResponseStatus(HttpStatus.SERVICE_UNAVAILABLE)
	public static class AsyncTimeoutException extends Exception {
	}
}