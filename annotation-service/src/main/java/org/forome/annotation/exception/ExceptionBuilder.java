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
import net.minidev.json.parser.ParseException;
import org.forome.annotation.struct.variant.Variant;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

public class ExceptionBuilder {

	private static final ExceptionFactory EXCEPTION_FACTORY = new ExceptionFactory();

	public static final String CODE_INVALID_CHROMOSOME = "invalid_chromosome";
	public static final String CODE_IO_ERROR = "io_error";

	public static AnnotatorException buildDatabaseException(DatabaseException cause) {
		return EXCEPTION_FACTORY.build("database_error", cause);
	}

	public static AnnotatorException buildExternalDatabaseException(Throwable cause) {
		return buildExternalDatabaseException(cause, null);
	}

	public static AnnotatorException buildExternalDatabaseException(String comment) {
		return buildExternalDatabaseException(null, comment);
	}

	public static AnnotatorException buildExternalDatabaseException(Throwable cause, String comment) {
		return EXCEPTION_FACTORY.build("external_database_error", comment, null, cause);
	}

	public static AnnotatorException buildExternalServiceException(Throwable cause) {
		return EXCEPTION_FACTORY.build("external_service_error", cause);
	}

	public static AnnotatorException buildExternalServiceException(Throwable cause, String serviceName, String comment) {
		return EXCEPTION_FACTORY.build("external_service_error", comment, new HashMap<String, Object>() {{
			put("service_name", serviceName);
		}}, cause);
	}

	public static AnnotatorException buildExternalServiceException(Throwable cause, String comment) {
		return EXCEPTION_FACTORY.build("external_service_error", comment, null, cause);
	}

	public static AnnotatorException buildNotMultipartRequestException() {
		return EXCEPTION_FACTORY.build("not_multipart_request");
	}

	public static AnnotatorException buildFileNotUploadException() {
		return EXCEPTION_FACTORY.build("file_not_upload");
	}

	public static AnnotatorException buildIOErrorException(IOException e) {
		return EXCEPTION_FACTORY.build(CODE_IO_ERROR, e);
	}

	public static AnnotatorException buildOperationException(Throwable cause) {
		return EXCEPTION_FACTORY.build("operation_error", cause);
	}

	public static AnnotatorException buildInvalidCredentialsException() {
		return EXCEPTION_FACTORY.build("invalid_credentials");
	}

	public static AnnotatorException buildInvalidVepJsonException(String value, Throwable cause) {
		return EXCEPTION_FACTORY.build("operation_error", null, new HashMap<String, Object>() {{
			put("value", value);
		}}, cause);
	}

	public static AnnotatorException buildLargeVcfFile(int maxSize) {
		return EXCEPTION_FACTORY.build("large_vcf_file", new HashMap<String, Object>() {{
			put("max_size", maxSize);
		}});
	}

	public static AnnotatorException buildInvalidVcfFile(Throwable cause) {
		return EXCEPTION_FACTORY.build("invalid_vcf_file", cause);
	}

	public static AnnotatorException buildNotEqualSamplesVcfAndFamFile(Collection<String> diff) {
		return EXCEPTION_FACTORY.build("not_equal_samples_vcf_and_fam_file", new HashMap<String, Object>() {{
			put("diff", diff);
		}});
	}

	public static AnnotatorException buildNotEqualSamplesVcfAndCnvFile() {
		return EXCEPTION_FACTORY.build("not_equal_samples_vcf_and_cnv_file");
	}

	public static AnnotatorException buildInvalidChromosome(String value) {
		return EXCEPTION_FACTORY.build(CODE_INVALID_CHROMOSOME, new HashMap<String, Object>() {{
			put("value", value);
		}});
	}

	public static AnnotatorException buildInvalidSequenceReference(Variant variant, String expected) {
		return EXCEPTION_FACTORY.build("invalid_sequence_reference", new HashMap<String, Object>() {{
			put("variant", variant.toString());
			put("expected", expected);
		}});
	}

	public static AnnotatorException buildNotUniqueValueException(String fieldName, Object fieldValue) {
		return EXCEPTION_FACTORY.build("not_unique_value", new HashMap<String, Object>() {{
			put("field_name", fieldName);
			put("field_value", fieldValue);
		}});
	}

	public static AnnotatorException buildInvalidValueException(String fieldName) {
		return EXCEPTION_FACTORY.build("invalid_value", new HashMap<String, Object>() {{
			put("field_name", fieldName);
		}});
	}

	public static AnnotatorException buildInvalidValueException(String fieldName, Object fieldValue) {
		return buildInvalidValueException(fieldName, fieldValue, null);
	}

	public static AnnotatorException buildInvalidValueException(String fieldName, Object fieldValue, String comment) {
		return EXCEPTION_FACTORY.build("invalid_value", comment, new HashMap<String, Object>() {{
			put("field_name", fieldName);
			put("field_value", fieldValue);
		}});
	}

	public static AnnotatorException buildInvalidValueJsonException(String fieldName, Throwable cause) {
		return EXCEPTION_FACTORY.build("invalid_value_json", null, new HashMap<String, Object>() {{
			put("field_name", fieldName);
		}}, cause);
	}

	public static AnnotatorException buildInvalidJsonException(ParseException cause) {
		return EXCEPTION_FACTORY.build("invalid_json", null, null, cause);
	}

	public static AnnotatorException buildInvalidInventoryException(Throwable cause) {
		return EXCEPTION_FACTORY.build("invalid_inventory", null, null, cause);
	}

	public static AnnotatorException buildInvalidValueInventoryException(String fieldName) {
		return buildInvalidValueInventoryException(fieldName, null);
	}

	public static AnnotatorException buildInvalidValueInventoryException(String fieldName, String comment) {
		return EXCEPTION_FACTORY.build("invalid_value_inventory", comment, new HashMap<String, Object>() {{
			put("field_name", fieldName);
		}}, null);
	}

	public static AnnotatorException buildInvalidOperation(String comment) {
		return EXCEPTION_FACTORY.build("invalid_operation", comment);
	}

	public static AnnotatorException buildServerTimeoutException() {
		return EXCEPTION_FACTORY.build("server_timeout");
	}

	public static AnnotatorException buildServerBusyException(String cause) {
		return EXCEPTION_FACTORY.build("server_busy", Collections.singletonMap("cause", cause));
	}

	public static AnnotatorException buildServerOverloadedException() {
		return EXCEPTION_FACTORY.build("server_overloaded");
	}

	public static AnnotatorException buildServerShutsDownException() {
		return EXCEPTION_FACTORY.build("server_shuts_down");
	}
}


