package ru.processtech.forome.annotation.network.component;

import com.infomaximum.database.domainobject.filter.HashFilter;
import ru.processtech.forome.annotation.database.entityobject.user.UserReadable;
import ru.processtech.forome.annotation.executionqueue.ExecutionTransaction;
import ru.processtech.forome.annotation.executionqueue.ReadableResource;
import ru.processtech.forome.annotation.executionqueue.ResourceProvider;

import java.util.Arrays;

public class AuthComponent {

	private final ReadableResource<UserReadable> userReadableResource;

	public AuthComponent(ResourceProvider resources) {
		userReadableResource = resources.getReadableResource(UserReadable.class);
	}

	public UserReadable auth(String login, String password, ExecutionTransaction transaction) {
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
