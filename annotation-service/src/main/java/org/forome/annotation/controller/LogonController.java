package org.forome.annotation.controller;

import com.infomaximum.querypool.Query;
import com.infomaximum.querypool.QueryTransaction;
import com.infomaximum.querypool.ResourceProvider;
import net.minidev.json.JSONObject;
import org.forome.annotation.Service;
import org.forome.annotation.controller.utils.ResponseBuilder;
import org.forome.annotation.database.entityobject.user.UserReadable;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.network.component.AuthComponent;
import org.forome.annotation.network.session.SessionService;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;

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

		return service.getQueryPool().execute(service.getDatabaseService().getDomainObjectSource(), new Query<String>() {

			private AuthComponent authComponent;

			@Override
			public void prepare(ResourceProvider resources) {
				authComponent = new AuthComponent(resources);
			}

			@Override
			public String execute(QueryTransaction transaction) {
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
