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

import org.forome.database.anotation.Entity;
import org.forome.database.anotation.Field;
import org.forome.database.anotation.HashIndex;
import org.forome.database.domainobject.DomainObject;

import javax.crypto.SecretKeyFactory;
import javax.crypto.spec.PBEKeySpec;
import java.security.NoSuchAlgorithmException;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.KeySpec;

@Entity(
		namespace = "org.forome.annotation",
		name = "User",
		fields = {
				@Field(name = "login", number = UserReadable.FIELD_LOGIN, type = String.class),
				@Field(name = "password_hash", number = UserReadable.FIELD_PASSWORD_HASH, type = byte[].class),
				@Field(name = "salt", number = UserReadable.FIELD_SALT, type = byte[].class),
		},
		hashIndexes = {
				@HashIndex(fields = { UserReadable.FIELD_LOGIN }),
		}
)
public class UserReadable extends DomainObject {

	public final static int FIELD_LOGIN = 0;
	public final static int FIELD_PASSWORD_HASH = 1;
	public final static int FIELD_SALT = 2;

	public UserReadable(long id) {
		super(id);
	}

	public String getLogin() {
		return getString(FIELD_LOGIN);
	}

	public byte[] getPasswordHash() {
		return get(FIELD_PASSWORD_HASH);
	}

	public byte[] getSalt() {
		return get(FIELD_SALT);
	}

	public static byte[] getSaltyPasswordHash(String password, byte[] salt) {
		if (password == null || salt == null) {
			return null;
		}
		try {
			SecretKeyFactory secretKeyFactory = SecretKeyFactory.getInstance("PBKDF2WithHmacSHA512");
			KeySpec ks = new PBEKeySpec(password.toLowerCase().toCharArray(), salt, 20000, 512);
			return secretKeyFactory.generateSecret(ks).getEncoded();
		} catch (InvalidKeySpecException | NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
	}
}
