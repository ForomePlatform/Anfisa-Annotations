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
import org.forome.querypool.QueryTransaction;
import org.forome.querypool.ReadableResource;
import org.forome.querypool.ResourceProvider;
import org.forome.annotation.database.entityobject.user.UserReadable;

import java.util.Arrays;

public class AuthComponent {

	private final ReadableResource<UserReadable> userReadableResource;

	public AuthComponent(ResourceProvider resources) {
		userReadableResource = resources.getReadableResource(UserReadable.class);
	}

	public UserReadable auth(String login, String password, QueryTransaction transaction) {
		UserReadable user = userReadableResource.find(new HashFilter(UserReadable.FIELD_LOGIN, login), transaction);
		if (user == null) return null;

		byte[] countSaltyPasswordHash = UserReadable.getSaltyPasswordHash(password, user.getSalt());
		if (Arrays.equals(countSaltyPasswordHash, user.getPasswordHash())) {
			return user;
		} else {
			return null;
		}
	}
}
