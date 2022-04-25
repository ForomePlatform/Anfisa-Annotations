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

package org.forome.annotation.annotator.main;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.concurrent.TimeUnit;

public class AnnotatorMainFork {

	public static void main(String[] args) throws IOException, InterruptedException, URISyntaxException {
		Path pathJava = Paths.get(System.getProperty("java.home")).resolve("bin").resolve("java").toAbsolutePath();
		if (!Files.exists(pathJava)) {
			throw new RuntimeException("java not found: " + pathJava);
		}

		String cmd = new StringBuilder()
				.append(pathJava.toString())
				.append(" -cp ").append(getPathJar())
				.append(" org.forome.annotation.annotator.main.AnnotatorMain ")
				.append(String.join(" ", args))
				.toString();

		System.out.println("Run process: " + cmd);

		ProcessBuilder processBuilder = new ProcessBuilder(cmd.split(" "));
		processBuilder.directory(Paths.get(".").toAbsolutePath().toFile());
		Process process = processBuilder.start();

		//Читаем stdout
		StringBuilder stdout = new StringBuilder();
		Thread threadStdOut = new Thread(() -> {
			try (InputStreamReader inputReader = new InputStreamReader(process.getInputStream())) {
				String line;
				try (BufferedReader bufferedReader = new BufferedReader(inputReader)) {
					while ((line = bufferedReader.readLine()) != null) {
						stdout.append(line).append(System.lineSeparator());
					}
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		});
//		threadStdOut.setDaemon(true);
		threadStdOut.start();

		//Читаем stderr
		StringBuilder stderr = new StringBuilder();
		Thread threadStdErr = new Thread(() -> {
			try (InputStreamReader inputReader = new InputStreamReader(process.getErrorStream())) {
				String line;
				try (BufferedReader bufferedReader = new BufferedReader(inputReader)) {
					while ((line = bufferedReader.readLine()) != null) {
						stderr.append(line).append(System.lineSeparator());
					}
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		});
//		threadStdErr.setDaemon(true);
		threadStdErr.start();


		//Немного ждем - убеждаемся, что все работам и выходим
		boolean complete = process.waitFor(7, TimeUnit.SECONDS);
		if (complete) {
			int exitCode = process.exitValue();
			if (exitCode == 0) {
				System.out.println("Process completed successfully");
			} else {
				System.out.println("The process ended with an error!!! ExitCode: " + exitCode);
				if (stdout.length() > 0) {
					System.out.println("stdout: ");
					System.out.println(stdout.toString());
				}
				if (stderr.length() > 0) {
					System.out.println("stderr: ");
					System.out.println(stderr.toString());
				}
			}
		}
	}

	private static Path getPathJar() throws URISyntaxException {
		URI jarLocationUri = AnnotatorMainFork.class.getProtectionDomain().getCodeSource().getLocation().toURI();
		Path pathJar = Paths.get(new File(jarLocationUri).getPath()).toAbsolutePath();
		if (!Files.exists(pathJar)) {
			throw new RuntimeException("jar file not found: " + pathJar);
		}
		return pathJar;
	}

}
