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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.Map;

public class RuntimeExec {

	public static class Result {

		public final int exitCode;
		public final String out;
		public final String outError;

		private Result(int exitCode, String out, String outError) {
			this.exitCode = exitCode;
			this.out = out;
			this.outError = outError;
		}
	}

	public static Result runCommand(String command) throws Exception{
		return runCommand(command, null);
	}

	public static Result runCommand(String command, Map<String, String> args) throws Exception{
		return runCommand(command, args, null, null);
	}

	public static Result runCommand(String command, Map<String, String> args, Map<String, String> environment, InputStream stdin) throws Exception{

		//Строим запускаемую команду
		String cmd = command;
		if (args!=null) {
			for(String key: args.keySet()) {
				cmd = cmd.replaceAll("\\$\\{" + key + "\\}", args.get(key));
			}
		}

		final ProcessBuilder processBuilder = new ProcessBuilder(cmd.split(" "));

		final Map<String, String> processBuilderEnvironment = processBuilder.environment();
		if (environment!=null) {
			for (String key: environment.keySet()) {
				processBuilderEnvironment.put(key, environment.get(key));
			}
		}

		//Выполняем
		Process process = processBuilder.start();

		//Возможно мы хотим писать во входной поток приложению
		if (stdin!=null) {
			OutputStream os = process.getOutputStream();
			int iByte;
			while ((iByte = stdin.read()) != -1) {
				os.write(iByte);
			}
			os.flush();
			os.close();
		}

		//Считываем выходной поток данных
		StringBuilder out = new StringBuilder();
		String line;
		BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
		while ((line = in.readLine()) != null) {
			out.append(line);
		}
		in.close();

		//Считываем выходной поток ошибочных данных
		StringBuilder outError = new StringBuilder();
		BufferedReader inError  = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		while ((line = inError.readLine()) != null) {
			outError.append(line);
		}
		inError.close();

		int exitCode = process.waitFor();//Ожидаем завершения
		return new Result(
				exitCode,
				out.toString(),
				outError.toString()
		);
	}
}
