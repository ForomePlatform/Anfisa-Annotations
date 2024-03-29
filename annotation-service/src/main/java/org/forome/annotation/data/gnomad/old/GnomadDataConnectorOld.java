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

package org.forome.annotation.data.gnomad.old;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableList;
import org.forome.annotation.data.DatabaseConnector;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.struct.SourceMetadata;
import org.forome.core.struct.Chromosome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Closeable;
import java.sql.*;
import java.util.*;
import java.util.stream.Collectors;

public class GnomadDataConnectorOld implements Closeable {

    private static final Logger log = LoggerFactory.getLogger(GnomadDataConnectorOld.class);

    public static final ImmutableList<String> ANCESTRIES = ImmutableList.of(
            "AFR",
            "AMR",
            "ASJ",
            "EAS",
            "FIN",
            "NFE",
            "OTH"
    );

    public static final ImmutableList<String> SEX = ImmutableList.of(
            "Male",
            "Female"
    );

    private static final ImmutableList<String> POP_GROUPS =
            new ImmutableList.Builder()
                    .addAll(ANCESTRIES)
                    .addAll(SEX)
                    .build();

    private static final ImmutableList<String> DATA_SUFFIXES =
            new ImmutableList.Builder()
                    .add("raw", "POPMAX")
                    .addAll(POP_GROUPS)
                    .build();

    private static final ImmutableList<String> KEY_COLUMNS = ImmutableList.of(
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT"
    );

    public static final ImmutableList<String> DATA_PREFIXES = ImmutableList.of(
            "AN",
            "AC",
            "Hom"
    );

    private static final ImmutableList<String> HOM = ANCESTRIES.stream().map(s -> "Hom_" + s).collect(ImmutableList.toImmutableList());

    private static final ImmutableList<String> AGGREGATE_DATA_COLUMNS = ImmutableList.of(
            "AN_Female + AN_Male as AN",
            "AC_Female + AC_Male as AC",
            String.join(" + ", HOM) + " as Hom",
            "AN_POPMAX",
            "AC_POPMAX",
            "POPMAX"
    );

    private static final ImmutableList<String> DATA_COLUMNS =
            DATA_PREFIXES.stream().flatMap(
                    a -> DATA_SUFFIXES.stream()
                            .filter(b -> !(a.equals("Hom") && b.equals("POPMAX")))
                            .map(b -> a + '_' + b)
            ).collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf));

    private static final String TABLE = "gnomad.VARIANTS";

    private static final ImmutableList<String> COLUMNS = new ImmutableList.Builder()
            .addAll(KEY_COLUMNS)
            .addAll(AGGREGATE_DATA_COLUMNS)
            .addAll(DATA_COLUMNS)
            .build();

    public class Result {

        public final Map<String, Object> columns;

        public Result(ResultSet resultSet) throws SQLException {
            Map<String, Object> modifiableColumns = new HashMap<>();
            ResultSetMetaData resultSetMetaData = resultSet.getMetaData();
            for (int i = 1; i <= resultSetMetaData.getColumnCount(); i++) {
                String columnName = resultSetMetaData.getColumnName(i);
                modifiableColumns.put(columnName, resultSet.getObject(i));
            }
            this.columns = Collections.unmodifiableMap(modifiableColumns);
        }
    }

    private final DatabaseConnector databaseConnector;

    public GnomadDataConnectorOld(DatabaseConnector databaseConnector) {
        this.databaseConnector = databaseConnector;
    }

    public List<SourceMetadata> getSourceMetadata(){
        return databaseConnector.getSourceMetadata();
    }

    public List<Result> getData(
            Chromosome chromosome,
            long position,
            String ref,
            String alt,
            String fromWhat,
            boolean exact) {

        String base_sql = String.format(
                "SELECT %s FROM %s WHERE CHROM = '%%s' and POS = %%s",
                String.join(", ", COLUMNS), TABLE
        );

        String sql = null;
        if (ref != null && alt != null) {
            if (exact) {
                sql = base_sql + " and REF = '%s' and ALT = '%s'";
            } else {
                sql = base_sql + " and (REF LIKE '%%%s%%' and ALT LIKE '%%%s%%')";
            }
        }

        if (fromWhat != null) {
            List<String> q = Arrays.asList(fromWhat.toLowerCase().split(","));
            String s;
            if (q.size() == 1) {
                if (q.get(0).contains("exome")) {
                    s = "e";
                } else if (q.get(0).contains("genome")) {
                    s = "g";
                } else {
                    s = q.get(0);
                }
            } else {
                throw new RuntimeException("Not support many s");
            }
            if (!(q.contains("e") && q.contains("g"))) {
                sql = String.format("%s and `SOURCE` = '%s'", sql, s);
            }
        }

        try (Connection connection = databaseConnector.createConnection()) {
            Statement statement = connection.createStatement();
            ResultSet resultSet = statement.executeQuery(String.format(sql, chromosome.getChar(), position, ref, alt));
            boolean isEmptyResultSet = !resultSet.next();

            if (!exact && isEmptyResultSet && ref != null && alt != null) {
                resultSet.close();
                statement.close();

                statement = connection.createStatement();
                resultSet = statement.executeQuery(
                        String.format(base_sql, chromosome.getChar(), position - 1)
                );
            }

            List<Result> results = new ArrayList<>();
            resultSet.beforeFirst();
            if (!exact) {
                while (resultSet.next()) {
                    String diff_ref_alt = diff(ref, alt);
                    if (Objects.equals(diff_ref_alt, diff(resultSet.getString(4), resultSet.getString(5)))
                            ||
                            diff3(resultSet.getString(4), resultSet.getString(5), diff_ref_alt)
                    ) {
                        results.add(new Result(resultSet));
                    }
                }
            } else {
                while (resultSet.next()) {
                    results.add(new Result(resultSet));
                }
            }

            resultSet.close();
            statement.close();

            return results;
        } catch (SQLException ex) {
            throw ExceptionBuilder.buildExternalDatabaseException(ex);
        }
    }

    @Override
    public void close()  {
        databaseConnector.close();
    }


    @VisibleForTesting
    public static String diff(String s1, String s2) {
        if (s1.equals(s2)) {
            return "";
        } else if (s2.contains(s1)) {
            int idx = s2.indexOf(s1);
            return s2.substring(0, idx) + s2.substring(idx + s1.length());
        } else if (s1.length() < s2.length()) {
            List<Character> x = new ArrayList<>();
            List<Character> y = new ArrayList<>();
            for (int i = 0; i < s2.length(); i++) {
                int j = x.size();
                if (j < s1.length() && s1.charAt(j) == s2.charAt(i)) {
                    x.add(s2.charAt(i));
                } else {
                    y.add(s2.charAt(i));
                }
            }
            if (!y.isEmpty()) {
                return y.stream().map(e -> e.toString()).collect(Collectors.joining());
            }
            return null;
        } else if (s2.length() < s1.length()) {
            return "-" + diff(s2, s1);
        } else {
            return null;
        }
    }

    public static boolean diff3(String s1, String s2, String d) {
        if (Strings.isNullOrEmpty(d)) {
            return s1.equals(s2);
        }
        if (d.charAt(0) == '-') {
            return diff3(s2, s1, d.substring(1));
        }
        if (s1.length() + d.length() != s2.length()) {
            return false;
        }
        for (int i=0; i<s1.length(); i++) {
            String x = s1.substring(0, i);
            String y = s1.substring(i);
            if ((x + d + y).equals(s2)) {
                return true;
            }
        }
        return false;
    }
}
