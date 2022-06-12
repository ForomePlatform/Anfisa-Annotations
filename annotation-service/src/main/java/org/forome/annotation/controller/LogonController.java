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

package org.forome.annotation.controller;

import org.forome.querypool.Query;
import org.forome.querypool.QueryTransaction;
import org.forome.querypool.ResourceProvider;
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
