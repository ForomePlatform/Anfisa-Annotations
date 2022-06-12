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
import org.forome.annotation.database.entityobject.user.UserEditable;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.network.component.UserEditableComponent;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;

import javax.servlet.http.HttpServletRequest;
import java.util.concurrent.CompletableFuture;

/**
 * http://localhost:8095/user/create?login=test&password=test66&session=...
 */
@Controller
@RequestMapping(value = { "/user", "/annotationservice/user" })
public class UserController {

	@RequestMapping(value = { "create" })
	public CompletableFuture<ResponseEntity> create(HttpServletRequest request) {
		Service service = Service.getInstance();

		String sessionId = request.getParameter("session");

		String login = request.getParameter("login");
		String password = request.getParameter("password");

		return service.getQueryPool().execute(service.getDatabaseService().getDomainObjectSource(), new Query<Long>() {

			private UserEditableComponent userEditableComponent;

			@Override
			public void prepare(ResourceProvider resources) {
				userEditableComponent = new UserEditableComponent(resources);
			}

			@Override
			public Long execute(QueryTransaction transaction) {
				Long requestUserId = service.getNetworkService().sessionService.checkSession(sessionId);
				if (requestUserId == null) {
					throw ExceptionBuilder.buildInvalidCredentialsException();
				}

				UserEditable user = userEditableComponent.create(login, password, transaction);
				return user.getId();
			}
		})
				.thenApply(userId -> {
					JSONObject out = new JSONObject();
					out.put("user_id", userId);
					return out;
				})
				.thenApply(out -> ResponseBuilder.build(out))
				.exceptionally(throwable -> ResponseBuilder.build(throwable));

	}
}
