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

package org.forome.annotation.utils;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class AppVersion {

	private static String _hashVersion = null;

	public static String getVersion() {
		if (_hashVersion == null) {
			String version = AppVersion.class.getPackage().getImplementationVersion();
			if (version == null) {
				//Возможно это запуск из проекта
				Path fileBuildGradle = Paths.get("build.gradle");
				if (!Files.exists(fileBuildGradle))  {
					String prjDir = System.getenv ("PRJ_DIR");
					if (prjDir != null) {
						fileBuildGradle = Paths.get (prjDir, "build.gradle");
					}
				}
				if (Files.exists(fileBuildGradle)) {
					try {
						try (Stream<String> stream = Files.lines(fileBuildGradle)) {
							String lineWithVersion = stream.filter(s -> s.trim().startsWith("version ")).findFirst().orElse(null);
							if (lineWithVersion != null) {
								Pattern pattern = Pattern.compile("'(.+?)'");
								Matcher matcher = pattern.matcher(lineWithVersion);
								matcher.find();
								version = matcher.group(1);
							}
						}
					} catch (Exception e) {
						throw new RuntimeException("Exception find version", e);
					}
				}
			}
			if (version == null) throw new RuntimeException("Not found version app");

			_hashVersion = version;
		}
		return _hashVersion;
	}

	public static String getVersionFormat() {
		String version = getVersion();
		String[] v = version.split("\\.");
		return new StringBuilder()
				.append(v[0]).append('.')
				.append(v[1]).append('.')
				.append(v[2])
				.toString();
	}

}
