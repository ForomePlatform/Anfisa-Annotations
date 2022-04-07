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

package org.forome.annotation.config.sshtunnel;

import net.minidev.json.JSONObject;

public class SshTunnelConfig {

	public final String host;
	public final int port;
	public final String user;
	public final String key;

	public SshTunnelConfig(JSONObject parse) {
		host = parse.getAsString("host");
		port = parse.getAsNumber("port").intValue();
		user = parse.getAsString("user");
		key = parse.getAsString("key");
	}

	@Override
	public String toString ()
	{
		return String.format ("%s@%s:%d", user, host, port);
	}
}
