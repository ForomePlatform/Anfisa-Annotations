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

package org.forome.annotation.network.component;

import org.forome.database.domainobject.filter.HashFilter;
import org.forome.querypool.EditableResource;
import org.forome.querypool.QueryTransaction;
import org.forome.querypool.ResourceProvider;
import org.forome.annotation.database.entityobject.user.UserEditable;
import org.forome.annotation.database.entityobject.user.UserReadable;
import org.forome.annotation.exception.ExceptionBuilder;

import java.util.regex.Pattern;

public class UserEditableComponent {

	private static final Pattern PATTERN_SECURITY = Pattern.compile("^(?=.*[0-9])(?=.*([a-z]|[A-Z]))(?=\\S+$).{6,}$");

	private final EditableResource<UserEditable> userEditableResource;

	public UserEditableComponent(ResourceProvider resources) {
		userEditableResource = resources.getEditableResource(UserEditable.class);
	}

	public UserEditable create(String login, String password, QueryTransaction transaction) {
		UserEditable user = userEditableResource.find(new HashFilter(UserReadable.FIELD_LOGIN, login), transaction);
		if (user != null) throw ExceptionBuilder.buildNotUniqueValueException("login", login);

		if (!PATTERN_SECURITY.matcher(password).matches()) {
			throw ExceptionBuilder.buildInvalidValueException("password");
		}

		user = userEditableResource.create(transaction);
		user.setLogin(login);
		user.setPassword(password);
		userEditableResource.save(user, transaction);

		return user;
	}
}
