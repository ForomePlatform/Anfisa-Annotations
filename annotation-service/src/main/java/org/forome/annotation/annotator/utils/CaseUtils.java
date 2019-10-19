package org.forome.annotation.annotator.utils;

import org.forome.annotation.struct.mcase.Cohort;
import org.forome.annotation.struct.mcase.MCase;
import org.forome.annotation.struct.mcase.Sample;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class CaseUtils {

    /**
     * sample["name"]      = sample_map.get(id, id)
     * sample['family']    = family
     * sample['id']        = id
     * sample['father']    = father
     * sample['mother']    = mother
     * sample['sex']       = int(sex)
     * sample['affected']  = (int(affected) == 2)
     */
    public static MCase parseFamFile(InputStream isFam, InputStream isFamSampleName, InputStream isCohorts) throws IOException {
        SortedMap<String, Sample> samples = new TreeMap<>();

        Map<String, String> sampleNameMap = new HashMap<>();
        if (isFamSampleName != null) {
            try (BufferedReader isBFamSampleName = new BufferedReader(new InputStreamReader(isFamSampleName))) {
                String line;
                while ((line = isBFamSampleName.readLine()) != null) {
                    line = line.trim();
                    if (line.length() == 0) continue;

                    String[] values = Arrays.stream(line.split("\t")).filter(s -> !s.trim().isEmpty()).toArray(String[]::new);
                    if (values.length != 4) {
                        throw new RuntimeException("Not support FamSampleName format");
                    }
                    String id = values[3];
                    String name = values[0];
                    sampleNameMap.put(id, name);
                }
            }
        }

        //Раскидываем маски семплов по когортам
        LinkedHashMap<Cohort, List<Pattern>> cohortMaskSamples = new LinkedHashMap<>();
        if (isCohorts != null) {
            try (BufferedReader isBFamCohorts = new BufferedReader(new InputStreamReader(isCohorts))) {
                String line;
                while ((line = isBFamCohorts.readLine()) != null) {
                    line = line.trim();
                    if (line.length() == 0) continue;

                    String[] values = Arrays.stream(line.split("\t")).filter(s -> !s.trim().isEmpty()).toArray(String[]::new);
                    if (values.length != 2) {
                        throw new RuntimeException("Not support file cohort format");
                    }

                    Cohort cohort = new Cohort(values[0]);
                    List<Pattern> maskSamples = Arrays.stream(values[1].split(",")).map(s -> s.trim())
                            .filter(s -> !s.isEmpty())
                            .map(s -> Pattern.compile(
                                    "^" + s.replaceAll("\\*", "(.*)") + "$"
                            ))
                            .collect(Collectors.toList());
                    cohortMaskSamples.put(cohort, maskSamples);
                }
            }
        }

        try (BufferedReader isBFam = new BufferedReader(new InputStreamReader(isFam))) {
            String line;
            while ((line = isBFam.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                String[] sl = line.split("\t");

                String id = sl[1];
                String name = sampleNameMap.getOrDefault(id, id);
                String family = sl[0];
                String father = sl[2];
                String mother = sl[3];
                int sex = Integer.parseInt(sl[4]);
                boolean affected = Integer.parseInt(sl[5]) == 2;

                Cohort cohort = null;
                for (Map.Entry<Cohort, List<Pattern>> entry : cohortMaskSamples.entrySet()) {
                    for (Pattern pattern: entry.getValue()) {
                        if (pattern.matcher(id).matches()) {
                            if (cohort!=null) {
                                throw new RuntimeException("Cohort: Not unique matches sample: " + id);
                            }

                            cohort = entry.getKey();
                            break;
                        }
                    }
                }

                Sample sample = new Sample(id, name, family, father, mother, sex, affected, cohort);
                samples.put(id, sample);
            }
        }

        return new MCase.Builder(samples, new ArrayList<>(cohortMaskSamples.keySet())).build();
    }

}
