package org.forome.annotation.network.component;

import com.infomaximum.database.domainobject.filter.HashFilter;
import com.infomaximum.querypool.EditableResource;
import com.infomaximum.querypool.QueryTransaction;
import com.infomaximum.querypool.ResourceProvider;
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
