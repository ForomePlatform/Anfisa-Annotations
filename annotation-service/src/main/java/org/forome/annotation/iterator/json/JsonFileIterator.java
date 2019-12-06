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

package org.forome.annotation.iterator.json;

import net.minidev.json.JSONObject;
import net.minidev.json.parser.JSONParser;
import net.minidev.json.parser.ParseException;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.struct.Chromosome;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class JsonFileIterator implements Iterator<JSONObject>, AutoCloseable {

	private final InputStream inputStream;
	private final BufferedReader bufferedReader;

	private JSONObject nextValue;

	public JsonFileIterator(Path pathVepJson) {
		this(getInputStream(pathVepJson), pathVepJson.getFileName().toString().endsWith(".gz"));
	}

	public JsonFileIterator(InputStream inputStream, boolean gzip) {
		this.inputStream = inputStream;
		if (gzip) {
			try {
				this.bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(inputStream)));
			} catch (IOException e) {
				throw ExceptionBuilder.buildIOErrorException(e);
			}
		} else {
			this.bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
		}

		nextValue = readNextValue();
	}

	@Override
	public boolean hasNext() {
		return (nextValue != null);
	}

	@Override
	public JSONObject next() {
		if (nextValue == null) {
			throw new NoSuchElementException();
		}

		JSONObject value = nextValue;
		while (true) {
			nextValue = readNextValue();
			if (nextValue != null
					&& Chromosome.CHR_M == Chromosome.of(nextValue.getAsString("seq_region_name"))
			) {
				continue;//Игнорируем митохондрии
			}
			break;
		}

		return value;
	}

	private JSONObject readNextValue() {
		String line = null;
		try {
			line = bufferedReader.readLine();
			if (line == null) {
				return null;
			} else {
				return (JSONObject) new JSONParser(JSONParser.DEFAULT_PERMISSIVE_MODE).parse(line);
			}
		} catch (IOException e) {
			throw ExceptionBuilder.buildIOErrorException(e);
		} catch (ParseException e) {
			throw ExceptionBuilder.buildInvalidVepJsonException(line, e);
		}
	}

	@Override
	public void close() throws IOException {
		bufferedReader.close();
		inputStream.close();
	}

	private static InputStream getInputStream(Path file) {
		try {
			return Files.newInputStream(file);
		} catch (IOException e) {
			throw ExceptionBuilder.buildIOErrorException(e);
		}
	}
}
