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

package org.forome.annotation.config.connector.base;

import net.minidev.json.JSONObject;

public abstract class DatabaseConfigConnector extends SshTunnelConfigConnector {

	public final String mysqlHost;
	public final int mysqlPort;
	public final String mysqlDatabase;
	public final String mysqlUser;
	public final String mysqlPassword;

	public DatabaseConfigConnector(JSONObject parse) {
		super(parse);

		JSONObject parseMysql = (JSONObject) parse.get("mysql");
		mysqlHost = parseMysql.getAsString("host");
		mysqlPort = parseMysql.getAsNumber("port").intValue();
		mysqlDatabase = parseMysql.getAsString("database");
		mysqlUser = parseMysql.getAsString("user");
		mysqlPassword = parseMysql.getAsString("password");
	}

	@Override
	public String toString ()
	{
		String host;
		if (sshTunnelConfig == null)
			host = mysqlHost;
		else
			host = super.toString () + "[" + mysqlHost + "]";
		String jdbcUrl = "jdbc:mysql://" +
						 host + ':' +
						 mysqlPort + '/' + mysqlDatabase +
						 "?user=" + mysqlUser +
						 "&password=" + mysqlPassword;
		return jdbcUrl;
	}
}
