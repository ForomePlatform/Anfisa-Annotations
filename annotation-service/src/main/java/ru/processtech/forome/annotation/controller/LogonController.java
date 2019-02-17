package ru.processtech.forome.annotation.controller;

import net.minidev.json.JSONObject;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;
import ru.processtech.forome.annotation.Service;
import ru.processtech.forome.annotation.controller.utils.ResponseBuilder;
import ru.processtech.forome.annotation.database.entityobject.user.UserReadable;
import ru.processtech.forome.annotation.exception.ExceptionBuilder;
import ru.processtech.forome.annotation.executionqueue.Execution;
import ru.processtech.forome.annotation.executionqueue.ExecutionTransaction;
import ru.processtech.forome.annotation.executionqueue.ResourceProvider;
import ru.processtech.forome.annotation.network.component.AuthComponent;
import ru.processtech.forome.annotation.network.session.SessionService;

import javax.servlet.http.HttpServletRequest;
import java.util.concurrent.CompletableFuture;

/**
 * http://localhost:8095/logon/login?login=admin&password=b82nfGl5sdg
 * http://34.218.20.85/annotationservice/logon/login?login=admin&password=b82nfGl5sdg
 */
@Controller
@RequestMapping(value = { "/logon", "/annotationservice/logon" })
public class LogonController {

	@RequestMapping(value = { "login" })
	public CompletableFuture<ResponseEntity> login(HttpServletRequest request) {
		Service service = Service.getInstance();

		String login = request.getParameter("login");
		String password = request.getParameter("password");

		return service.getExecutionQueue().execute(service, new Execution<String>() {

			private AuthComponent authComponent;

			@Override
			public void prepare(ResourceProvider resources) {
				authComponent = new AuthComponent(resources);
			}

			@Override
			public String execute(ExecutionTransaction transaction) {
				UserReadable userReadable = authComponent.auth(login, password, transaction);
				if (userReadable == null) throw ExceptionBuilder.buildInvalidCredentialsException();

				return service.getNetworkService().sessionService.buildSession(userReadable);
			}
		}).thenApply(sessionId -> {
			JSONObject out = new JSONObject();
			out.put("session", sessionId);
			out.put("expire", System.currentTimeMillis() + SessionService.SESSION_TIMEOUT.toMillis());
			return out;
		})
				.thenApply(out -> ResponseBuilder.build(out))
				.exceptionally(throwable -> ResponseBuilder.build(throwable));

	}
}
