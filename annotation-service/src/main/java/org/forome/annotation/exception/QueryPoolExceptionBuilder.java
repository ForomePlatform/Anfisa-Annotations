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

package org.forome.annotation.exception;

import org.forome.database.exception.DatabaseException;
import org.forome.querypool.ExceptionBuilder;

public class QueryPoolExceptionBuilder implements ExceptionBuilder {

	@Override
	public RuntimeException buildDatabaseException(DatabaseException e) {
		return org.forome.annotation.exception.ExceptionBuilder.buildDatabaseException(e);
	}

	@Override
	public RuntimeException buildServerBusyException(String cause) {
		return org.forome.annotation.exception.ExceptionBuilder.buildServerBusyException(cause);
	}

	@Override
	public RuntimeException buildServerOverloadedException() {
		return org.forome.annotation.exception.ExceptionBuilder.buildServerOverloadedException();
	}

	@Override
	public RuntimeException buildServerShutsDownException() {
		return org.forome.annotation.exception.ExceptionBuilder.buildServerShutsDownException();
	}
}
