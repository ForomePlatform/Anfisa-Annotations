package org.forome.annotation.connector.ref;

import org.forome.annotation.connector.DatabaseConnector;
import org.forome.annotation.exception.ExceptionBuilder;
import org.forome.annotation.struct.Chromosome;
import org.forome.annotation.utils.NucleotideUtils;

import java.io.Closeable;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class RefDataConnector implements Closeable {

    private static final long GENE_BUCKET_SIZE = 1000000L;

    private final DatabaseConnector databaseConnector;

    public RefDataConnector(DatabaseConnector databaseConnector) {
        this.databaseConnector = databaseConnector;
    }

    public String getRef(Chromosome chromosome, int start, int end) {
        String sql = String.format(
                "SELECT Ref, Pos FROM util.hg19 WHERE Chrom = %s AND Pos between %s and %s ORDER BY Pos ASC",
                chromosome.getChar(), Math.min(start, end), Math.max(start, end)
        );

        StringBuilder ref = new StringBuilder();
        int lastPos = start - 1;
        try (Connection connection = databaseConnector.createConnection()) {
            try (Statement statement = connection.createStatement()) {
                try (ResultSet resultSet = statement.executeQuery(sql)) {
                    while (resultSet.next()) {
                        int pos = resultSet.getInt("Pos");
                        //Валидируем порядок
                        if (pos != ++lastPos) {
                            throw new RuntimeException("Bad sequence, sql: " + sql);
                        }
                        String sNucleotide = resultSet.getString("Ref");
                        if (sNucleotide.length() != 1) {
                            throw new RuntimeException("Bad ref, sql: " + sql + ", ref: " + sNucleotide);
                        }
                        char nucleotide = sNucleotide.toUpperCase().charAt(0);
                        if (!NucleotideUtils.validation(nucleotide)) {
                            throw new RuntimeException("Not valid nucleotide, sql: " + sql + ", nucleotide: " + nucleotide);
                        }

                        ref.append(nucleotide);
                    }
                }
            }
        } catch (SQLException ex) {
            throw ExceptionBuilder.buildExternalDatabaseException(ex);
        }
        if (ref.length() != Math.max(start, end) - Math.min(start, end) + 1) {
            throw new RuntimeException("Bad len ref, sql: " + sql + ", ref: " + ref);
        }
        return ref.toString();
    }

    @Override
    public void close() {
        databaseConnector.close();
    }
}
