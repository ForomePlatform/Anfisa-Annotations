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

package org.forome.annotation.database.entityobject.user;

import org.forome.database.domainobject.DomainObjectEditable;
import org.forome.annotation.utils.RandomUtils;

public class UserEditable extends UserReadable implements DomainObjectEditable {

	public UserEditable(long id) {
		super(id);
	}

	public void setLogin(String login) {
		set(FIELD_LOGIN, login);
	}

	public void setPassword(String password) {
		byte[] salt = null;
		if (password != null) {
			salt = new byte[16];
			RandomUtils.SECURE_RANDOM.nextBytes(salt);
		}
		set(FIELD_SALT, salt);
		set(FIELD_PASSWORD_HASH, getSaltyPasswordHash(password, salt));
	}

}
