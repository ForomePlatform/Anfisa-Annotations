package org.forome.annotation.connector.anfisa;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import org.forome.annotation.annotator.utils.SampleUtils;
import org.forome.annotation.connector.anfisa.struct.*;
import org.forome.annotation.connector.beacon.BeaconConnector;
import org.forome.annotation.connector.clinvar.ClinvarConnector;
import org.forome.annotation.connector.clinvar.struct.ClinvarResult;
import org.forome.annotation.connector.clinvar.struct.ClinvarVariantSummary;
import org.forome.annotation.connector.conservation.ConservationConnector;
import org.forome.annotation.connector.conservation.struct.Conservation;
import org.forome.annotation.connector.gnomad.GnomadConnector;
import org.forome.annotation.connector.gnomad.struct.GnomadResult;
import org.forome.annotation.connector.gtf.GTFConnector;
import org.forome.annotation.connector.gtf.struct.GTFRegion;
import org.forome.annotation.connector.hgmd.HgmdConnector;
import org.forome.annotation.connector.liftover.LiftoverConnector;
import org.forome.annotation.connector.spliceai.SpliceAIConnector;
import org.forome.annotation.connector.spliceai.struct.SpliceAIResult;
import org.forome.annotation.exception.AnnotatorException;
import org.forome.annotation.struct.Chromosome;
import org.forome.annotation.struct.Position;
import org.forome.annotation.struct.Sample;
import org.forome.annotation.struct.variant.Variant;
import org.forome.annotation.struct.variant.VariantVCF;
import org.forome.annotation.struct.variant.cnv.VariantCNV;
import org.forome.annotation.utils.AppVersion;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AnfisaConnector implements AutoCloseable {

    private final static Logger log = LoggerFactory.getLogger(AnfisaConnector.class);

    private static final Map<String, String> trustedSubmitters = new HashMap<String, String>() {{
        put("lmm", "Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine");
        put("gene_dx", "GeneDx");
    }};

    private static final Map<String, String> proteins_3_to_1 = new HashMap<String, String>() {{
        put("Ala", "A");
        put("Arg", "R");
        put("Asn", "N");
        put("Asp", "D");
        put("Cys", "C");
        put("Gln", "Q");
        put("Glu", "E");
        put("Gly", "G");
        put("His", "H");
        put("Ile", "I");
        put("Leu", "L");
        put("Lys", "K");
        put("Met", "M");
        put("Phe", "F");
        put("Pro", "P");
        put("Ser", "S");
        put("Thr", "T");
        put("Trp", "W");
        put("Tyr", "Y");
        put("Val", "V");

    }};

    public final GnomadConnector gnomadConnector;
    public final SpliceAIConnector spliceAIConnector;
    public final ConservationConnector conservationConnector;
    public final HgmdConnector hgmdConnector;
    public final ClinvarConnector clinvarConnector;
    private final LiftoverConnector liftoverConnector;
    private final GTFConnector gtfConnector;

    public AnfisaConnector(
            GnomadConnector gnomadConnector,
            SpliceAIConnector spliceAIConnector,
            ConservationConnector conservationConnector,
            HgmdConnector hgmdConnector,
            ClinvarConnector clinvarConnector,
            LiftoverConnector liftoverConnector,
            GTFConnector gtfConnector,
            Thread.UncaughtExceptionHandler uncaughtExceptionHandler
    ) throws IOException {
        this.gnomadConnector = gnomadConnector;
        this.spliceAIConnector = spliceAIConnector;
        this.conservationConnector = conservationConnector;
        this.hgmdConnector = hgmdConnector;
        this.clinvarConnector = clinvarConnector;
        this.liftoverConnector = liftoverConnector;

        this.gtfConnector = gtfConnector;
    }

    public AnfisaResult build(
            String caseSequence,
            AnfisaInput anfisaInput,
            Variant variant, JSONObject vepJson
    ) {
        Record record = new Record();
        AnfisaExecuteContext context = new AnfisaExecuteContext(
                anfisaInput, variant, vepJson
        );

        AnfisaResultFilters filters = new AnfisaResultFilters();
        AnfisaResultData data = new AnfisaResultData();
        AnfisaResultView view = new AnfisaResultView();

        data.version = AppVersion.getVersionFormat();

        callGnomAD(context, variant, anfisaInput.samples, vepJson, filters);
        callSpliceai(data, filters, variant, anfisaInput.samples, vepJson);
        callHgmd(record, context, filters, data);
        callClinvar(context, record, variant.chromosome.getChar(), variant.start, variant.end, anfisaInput.samples, filters, data, view, vepJson);
        callBeacon(variant, anfisaInput.samples, vepJson, data);
        GtfAnfisaResult gtfAnfisaResult = callGtf(variant.start, variant.end, vepJson);
        callQuality(filters, variant, anfisaInput.samples);

        filters.severity = getSeverity(vepJson);
        filters.alts = alt_list(variant, anfisaInput.samples, vepJson);

        String proband = SampleUtils.getProband(anfisaInput.samples);
        if (proband != null) {
            String mother = anfisaInput.samples.get(proband).mother;
            String father = anfisaInput.samples.get(proband).father;
            data.zygosity = new HashMap<>();
            filters.altZygosity = new HashMap<>();
            for (Map.Entry<String, Sample> entry : anfisaInput.samples.entrySet()) {
                String name = entry.getValue().name;
                int sex = entry.getValue().sex;
                String label;
                if (entry.getKey().equals(proband)) {
                    label = String.format("proband [%s]", name);
                } else if (entry.getKey().equals(mother)) {
                    label = String.format("mother [%s]", name);
                } else if (entry.getKey().equals(father)) {
                    label = String.format("father [%s]", name);
                } else {
                    label = entry.getKey();
                }

                Integer zyg = sampleHasVariant(vepJson, variant, anfisaInput.samples, entry.getValue());
                data.zygosity.put(entry.getKey(), zyg);
                Integer modified_zygosity = (!variant.chromosome.getChar().equals("X") || sex == 2 || (zyg != null && zyg == 0)) ? zyg : (Integer) 2;
                filters.altZygosity.put(entry.getKey(), modified_zygosity);
                if (zyg != null && zyg > 0) {
                    filters.has_variant.add(label);
                }
            }
        }

        List<Object> d = getDistanceFromExon(gtfAnfisaResult, vepJson, "worst");
        filters.distFromExon = d.stream()
                .filter(o -> (o instanceof Number))
                .map(o -> ((Number) o).longValue())
                .min(Long::compareTo).orElse(0L);

        filters.chromosome = variant.chromosome.getChromosome();

        data.assemblyName = vepJson.getAsString("assembly_name");
        data.end = vepJson.getAsNumber("end").longValue();
        data.regulatoryFeatureConsequences = (JSONArray) vepJson.get("regulatory_feature_consequences");
        data.motifFeatureConsequences = (JSONArray) vepJson.get("motif_feature_consequences");
        data.intergenicConsequences = (JSONArray) vepJson.get("intergenic_consequences");
        data.start = vepJson.getAsNumber("start").longValue();
        data.mostSevereConsequence = vepJson.getAsString("most_severe_consequence");
        data.alleleString = vepJson.getAsString("allele_string");
        data.seqRegionName = vepJson.getAsString("seq_region_name");
        data.colocatedVariants = (JSONArray) vepJson.get("colocated_variants");
        data.input = vepJson.getAsString("input");
        data.label = getLabel(context);
        data.transcriptConsequences = (JSONArray) vepJson.get("transcript_consequences");
        data.id = vepJson.getAsString("id");
        data.strand = (vepJson.containsKey("strand")) ? vepJson.getAsNumber("strand").longValue() : null;
        data.variantClass = (vepJson.containsKey("variant_class")) ? vepJson.getAsString("variant_class") : null;

        data.colorCode = getColorCode(vepJson, data, record, filters);

        data.distFromBoundaryCanonical = gtfAnfisaResult.distFromBoundaryCanonical;
        data.regionCanonical = gtfAnfisaResult.regionCanonical;
        data.distFromBoundaryWorst = gtfAnfisaResult.distFromBoundaryWorst;
        data.regionWorst = gtfAnfisaResult.regionWorst;

        if (variant instanceof VariantCNV) {
            VariantCNV variantCNV = (VariantCNV) variant;
            data.cnvGT = variantCNV.genotypes.values().stream()
                    .collect(Collectors.toMap(item -> item.sampleName, item -> item.gt));
            filters.cnvLO = view.bioinformatics.cnvLO =
                    variantCNV.getGenotype(SampleUtils.getProband(variantCNV.genotypes.keySet())).lo;
        }

        createGeneralTab(context, data, filters, view, variant.start, variant.end, vepJson, caseSequence, variant, anfisaInput.samples);
        createQualityTab(view, variant, anfisaInput.samples);
        createGnomadTab(context, variant.chromosome.getChar(), variant, anfisaInput.samples, vepJson, filters, data, view);
        createDatabasesTab(vepJson, record, data, view);
        createPredictionsTab(vepJson, view);
        createBioinformaticsTab(gtfAnfisaResult, context, data, view);

        return new AnfisaResult(filters, data, view);
    }

    private static Integer sampleHasVariant(JSONObject json, Variant variant, Map<String, Sample> samples, Sample sample) {
        if (variant == null || !(variant instanceof VariantVCF)) {
            return 0;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        Integer idx = null;
        if (sample.name.toLowerCase().equals("proband")) {
            idx = 0;
        } else if (sample.name.toLowerCase().equals("mother")) {
            idx = 1;
        } else if (sample.name.toLowerCase().equals("father")) {
            idx = 2;
        }

        Genotype oGenotype;
        if (idx == null) {
            oGenotype = variantContext.getGenotype(sample.id);
        } else {
            /**
             genotypes = self.get_genotypes()
             if (len(genotypes) <= idx):
             return False
             genotype = genotypes[idx]
             */
            throw new RuntimeException("Not implemented");
        }

        String genotype;
        if (oGenotype.isCalled()) {
            genotype = oGenotype.getGenotypeString();
        } else if (oGenotype.getType() == GenotypeType.UNAVAILABLE) {
            //Не имеет альтернативного аллеля
            return 0;
        } else if (oGenotype.getType() == GenotypeType.NO_CALL) {
            //Генотип не может быть определен из-за плохого качества секвенирования
            return null;
        } else {
            throw new RuntimeException("Unknown state");
        }

        Set<String> set1 = Arrays.stream(genotype.split("/")).collect(Collectors.toSet());
        Set<String> set2 = new HashSet<>(alt_list(variant, samples, json));
        if (!Collections.disjoint(set1, set2)) {
            return 3 - set1.size();
        }
        return 0;
    }

    private void callHgmd(Record record, AnfisaExecuteContext anfisaExecuteContext, AnfisaResultFilters filters, AnfisaResultData data) {
        Variant variant = anfisaExecuteContext.variant;
        List<String> accNums = hgmdConnector.getAccNum(variant.chromosome.getChar(), variant.start, variant.end);
        if (accNums.size() > 0) {
            HgmdConnector.Data hgmdData = hgmdConnector.getDataForAccessionNumbers(accNums);
            record.hgmdData = hgmdData;

            data.hgmd = String.join(",", accNums);
            List<Long[]> hg38 = hgmdConnector.getHg38(accNums);

            data.hgmdHg38 = hg38.stream().map(longs -> String.format("%s-%s", longs[0], longs[1])).collect(Collectors.joining(", "));
            List<String> tags = hgmdData.hgmdPmidRows.stream().map(hgmdPmidRow -> hgmdPmidRow.tag).collect(Collectors.toList());
            filters.hgmdBenign = (tags.size() == 0);
        }
    }

    private void callClinvar(AnfisaExecuteContext context, Record record, String _chromosome, long start, long end, Map<String, Sample> samples, AnfisaResultFilters filters, AnfisaResultData data, AnfisaResultView view, JSONObject json) {
        Variant variant = context.variant;
        Chromosome chromosome = variant.chromosome;

        List<ClinvarResult> clinvarResults;
        if (isSnv(json)) {
            clinvarResults = clinvarConnector.getData(_chromosome, start, end, alt_list(variant, samples, json));
        } else {
            clinvarResults = clinvarConnector.getExpandedData(_chromosome, start);
        }
        record.clinvarResults = clinvarResults;
        if (!clinvarResults.isEmpty()) {

            String[] variants = clinvarResults.stream().map(clinvarResult -> {
                return String.format("%s %s>%s",
                        vstr(chromosome.getChromosome(), new Position(clinvarResult.start, clinvarResult.end)),
                        clinvarResult.referenceAllele, clinvarResult.alternateAllele
                );
            }).toArray(String[]::new);

            List<String> significance = new ArrayList<>();
            Map<String, String> submissions = new HashMap<>();
            for (ClinvarResult clinvarResult : clinvarResults) {
                significance.addAll(Arrays.asList(clinvarResult.clinicalSignificance.split("/")));
                submissions.putAll(clinvarResult.submitters);
            }

            List<String> idList = clinvarResults.stream().flatMap(it -> {
                return Lists.newArrayList(it.phenotypeIDs, it.otherIDs).stream();
            }).collect(Collectors.toList());
            for (String id : idList) {
                if (id.indexOf(":") != -1) {
                    continue;
                }
                //TODO not implemented
            }

            data.clinVar = clinvarResults.stream()
                    .map(clinvarResult -> clinvarResult.variationID)
                    .map(it -> Long.parseLong(it))
                    .toArray(Long[]::new);
            data.clinvarSubmitters = new HashMap<String, String>() {{
                for (ClinvarResult clinvarResult : clinvarResults) {
                    for (Map.Entry<String, String> entry : clinvarResult.submitters.entrySet()) {
                        put(encodeToAscii(entry.getKey()), entry.getValue());
                    }
                }
            }};

            view.databases.clinVar = clinvarResults.stream()
                    .map(clinvarResult -> clinvarResult.variationID)
                    .map(it -> String.format("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/", it))
                    .toArray(String[]::new);
            data.clinvarVariants = variants;
            view.databases.clinVarSubmitters = data.clinvarSubmitters.entrySet().stream().map(entry -> {
                return String.format("%s: %s", encodeToAscii(entry.getKey()), entry.getValue());
            }).sorted().toArray(String[]::new);
            data.clinvarSignificance = significance.toArray(new String[significance.size()]);
            data.clinvarPhenotypes = clinvarResults.stream()
                    .map(clinvarResult -> clinvarResult.phenotypeList)
                    .toArray(String[]::new);

            filters.clinvarBenign = (significance.stream().filter(s -> (s.toLowerCase().indexOf("benign") == -1)).count() == 0);

            Boolean benign = null;
            for (String submitter : trustedSubmitters.keySet()) {
                String fullName = trustedSubmitters.get(submitter);
                if (submissions.containsKey(fullName)) {
                    String prediction = submissions.get(fullName).toLowerCase();
                    data.setField(submitter, prediction);
                    if (!prediction.contains("benign")) {
                        benign = false;
                    } else if (benign == null) {
                        benign = true;
                    }
                }
            }
            filters.clinvarTrustedBenign = Optional.ofNullable(benign);
        }

        ClinvarVariantSummary clinvarVariantSummary = clinvarConnector.getDataVariantSummary(chromosome, start, end);
        if (clinvarVariantSummary != null) {
            view.databases.clinvarReviewStatus = clinvarVariantSummary.reviewStatus.text;
            filters.clinvarReviewStatus = clinvarVariantSummary.reviewStatus;
            view.databases.numClinvarSubmitters = clinvarVariantSummary.numberSubmitters;
            filters.numClinvarSubmitters = clinvarVariantSummary.numberSubmitters;
            view.databases.clinvarAcmgGuidelines = filters.clinvarAcmgGuidelines =
                    clinvarVariantSummary.guidelineTypes.stream()
                            .map(guideline -> guideline.getStrValue()).toArray(String[]::new);
        }
    }

    private void callBeacon(Variant variant, Map<String, Sample> samples, JSONObject json, AnfisaResultData data) {
        List<String> alts = alt_list(variant, samples, json);
        data.beaconUrls = alts.stream()
                .map(alt ->
                        BeaconConnector.getUrl(
                                variant.chromosome.getChar(),
                                json.getAsNumber("start").longValue(),
                                getRef(variant, json), alt
                        )
                )
                .toArray(String[]::new);
    }

    private GtfAnfisaResult callGtf(long start, long end, JSONObject json) {
        GtfAnfisaResult gtfAnfisaResult = new GtfAnfisaResult();

        //TODO Ulitin V. Отличие от python-реализации
        //Дело в том, что в оригинальной версии используется set для позиции, но в коде ниже используется итерация этому
        //списку и в конечном итоге это вляет на значение поля region - судя по всему это потенциальный баг и
        //необходима консультация с Михаилом
        List<Long> pos = new ArrayList<>();
        pos.add(start);
        if (start != end) {
            pos.add(end);
        }

        List<String> transcriptKinds = Lists.newArrayList("canonical", "worst");
        for (String kind : transcriptKinds) {
            List<String> transcripts;
            if ("canonical".equals(kind)) {
                transcripts = getCanonicalTranscripts(json).stream()
                        .filter(jsonObject -> "Ensembl".equals(jsonObject.getAsString("source")))
                        .map(jsonObject -> jsonObject.getAsString("transcript_id"))
                        .collect(Collectors.toList());
            } else if ("worst".equals(kind)) {
                transcripts = getMostSevereTranscripts(json).stream()
                        .filter(jsonObject -> "Ensembl".equals(jsonObject.getAsString("source")))
                        .map(jsonObject -> jsonObject.getAsString("transcript_id"))
                        .collect(Collectors.toList());
            } else {
                throw new RuntimeException("Unknown Transcript Kind: " + kind);
            }
            if (transcripts.isEmpty()) {
                continue;
            }

            List<Object[]> distances = new ArrayList<>();
            Long dist = null;
            String region = null;
            Long index = null;
            Integer n = null;
            for (String t : transcripts) {
                dist = null;
                for (Long p : pos) {
                    Object[] result = gtfConnector.lookup(p, t);
                    if (result == null) {
                        continue;
                    }
                    long d = (long) result[0];

                    GTFRegion gtfRegion = (GTFRegion) result[1];
                    region = gtfRegion.region;
                    if (gtfRegion.indexRegion != null) {
                        index = gtfRegion.indexRegion;
                        n = (int) result[2];
                    }

                    if (dist == null || d < dist) {
                        dist = d;
                    }
                }
                distances.add(
                        new Object[]{
                                dist, region, index, n
                        }
                );
            }

            if ("canonical".equals(kind)) {
                gtfAnfisaResult.distFromBoundaryCanonical = distances;
                gtfAnfisaResult.regionCanonical = region;
            } else if ("worst".equals(kind)) {
                gtfAnfisaResult.distFromBoundaryWorst = distances;
                gtfAnfisaResult.regionWorst = region;
            } else {
                throw new RuntimeException("Unknown Transcript Kind: " + kind);
            }
        }
        return gtfAnfisaResult;
    }

    private static void callQuality(AnfisaResultFilters filters, Variant variant, Map<String, Sample> samples) {
        filters.minGq = getMinGQ(variant, samples);
        filters.probandGq = getProbandGQ(variant, samples);
        if (variant != null && variant instanceof VariantVCF) {
            VariantContext variantContext = ((VariantVCF) variant).variantContext;

            CommonInfo commonInfo = variantContext.getCommonInfo();
            filters.qd = toPrimitiveDouble(commonInfo.getAttribute("QD"));
            filters.fs = toPrimitiveDouble(commonInfo.getAttribute("FS"));
            filters.mq = toDouble(commonInfo.getAttribute("MQ"));

            if (variantContext.isFiltered()) {
                filters.filters = new ArrayList<>(variantContext.getFilters());
            } else {
                filters.filters = Lists.newArrayList("PASS");
            }
        }
    }


    private void callGnomAD(AnfisaExecuteContext context, Variant variant, Map<String, Sample> samples, JSONObject response, AnfisaResultFilters filters) {
        Double af = null;
        Double _af = null;
        Double emAfPb = null;
        Double gmAfPb = null;
        Double _afPb = null;

        GnomadResult.Popmax popmax = null;
        GnomadResult.Popmax widePopmax = null;

        for (String alt : alt_list(variant, samples, response)) {
            GnomadResult gnomadResult = getGnomadResult(variant, response, alt);
            if (gnomadResult == null) {
                continue;
            }
            if (gnomadResult.exomes != null) {
                af = gnomadResult.exomes.af;
                if (isProbandHasAllele(variant, samples, alt)) {
                    emAfPb = Math.min((emAfPb != null) ? emAfPb : af, af);
                }
            }
            if (gnomadResult.genomes != null) {
                af = gnomadResult.genomes.af;
                if (isProbandHasAllele(variant, samples, alt)) {
                    gmAfPb = Math.min((gmAfPb != null) ? gmAfPb : af, af);
                }
            }

            af = gnomadResult.overall.af;
            if (isProbandHasAllele(variant, samples, alt)) {
                _afPb = Math.min((_afPb != null) ? _afPb : af, af);
            }

            if (_af == null || af < _af) {
                _af = af;
                popmax = gnomadResult.popmax;
                widePopmax = gnomadResult.widePopmax;
            }
        }

        context.gnomadAfFam = _af;
        filters.gnomadAfPb = _afPb;

        filters.gnomadPopmax = popmax;
        filters.gnomadWidePopmax = widePopmax;
    }

    private GnomadResult getGnomadResult(Variant variant, JSONObject response, String alt) {
        try {
            return gnomadConnector.request(
                    variant.chromosome.getChar(),
                    Math.min(response.getAsNumber("start").longValue(), response.getAsNumber("end").longValue()),
                    getRef(variant, response), alt
            ).get();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof AnnotatorException) {
                throw (AnnotatorException) cause;
            } else {
                throw new RuntimeException(cause);
            }
        }
    }

    private void callSpliceai(AnfisaResultData data, AnfisaResultFilters filters, Variant variant, Map<String, Sample> samples, JSONObject json) {
        SpliceAIResult spliceAIResult = spliceAIConnector.getAll(
                variant.chromosome.getChar(),
                lowest_coord(json),
                ref(json, variant),
                alt_list(variant, samples, json)
        );
        data.spliceAI = spliceAIResult.dict_sql;
        filters.spliceAltering = spliceAIResult.cases;
        filters.spliceAiDsmax = spliceAIResult.max_ds;
    }

    private void createGeneralTab(AnfisaExecuteContext context, AnfisaResultData data, AnfisaResultFilters filters, AnfisaResultView view, long start, long end, JSONObject json, String caseSequence, Variant variant, Map<String, Sample> samples) {
        view.general.genes = getGenes(json).stream().toArray(String[]::new);
        view.general.hg19 = str(context);
        view.general.hg38 = getStrHg38Coordinates(context);

        if (isSnv(json)) {
            data.ref = getRef(variant, json);
            data.alt = altString(variant, samples, json);
        } else {
            view.general.ref = getRef(variant, json);
            view.general.alt = altString(variant, samples, json);
        }

        List<String>[] cPosTpl = getPosTpl(json, "c");
        view.general.cposWorst = cPosTpl[0];
        view.general.cposCanonical = cPosTpl[1];
        view.general.cposOther = cPosTpl[2];

        List<String>[] pPosTpl = getPosTpl(json, "p");
        view.general.pposWorst = pPosTpl[0];
        view.general.pposCanonical = pPosTpl[1];
        view.general.pposOther = pPosTpl[2];

        Object[] gGenotypes = getGenotypes(variant, samples);
        view.general.probandGenotype = (String) gGenotypes[0];
        view.general.maternalGenotype = (String) gGenotypes[1];
        view.general.paternalGenotype = (String) gGenotypes[2];

        view.general.worstAnnotation = data.mostSevereConsequence;
        List<String> consequenceTerms = getFromCanonicalTranscript(json, "consequence_terms");
        String canonicalAnnotation = getMostSevere(consequenceTerms);
        if (consequenceTerms.size() > 1) {
            String finalCanonicalAnnotation = canonicalAnnotation;
            List<String> otherTerms = consequenceTerms.stream()
                    .filter(s -> !s.equals(finalCanonicalAnnotation))
                    .collect(Collectors.toList());
            canonicalAnnotation = String.format("%s [%s]", canonicalAnnotation, String.join(", ", otherTerms));
        }
        view.general.canonicalAnnotation = canonicalAnnotation;

        view.general.spliceRegion = getFromTranscripts(json, "spliceregion", "all");
        view.general.geneSplicer = getFromTranscripts(json, "genesplicer", "all");

        List<JSONObject> transcripts = getMostSevereTranscripts(json);
        view.general.refseqTranscriptWorst = getFromTranscripts(transcripts, "transcript_id", "RefSeq");
        view.general.ensemblTranscriptsWorst = getFromTranscripts(transcripts, "transcript_id", "Ensembl");

        transcripts = getCanonicalTranscripts(json);
        view.general.refseqTranscriptCanonical = getFromTranscripts(transcripts, "transcript_id", "RefSeq");
        view.general.ensemblTranscriptsCanonical = getFromTranscripts(transcripts, "transcript_id", "Ensembl");

        view.general.variantExonWorst = getFromWorstTranscript(json, "exon");
        view.general.variantIntronWorst = getFromWorstTranscript(json, "intron");
        view.general.variantExonCanonical = getFromCanonicalTranscript(json, "exon");
        view.general.variantIntronCanonical = getFromCanonicalTranscript(json, "intron");

        String[] intronOrExonCanonical = getIntronOrExon(json, "canonical");
        data.variantExonIntronCanonical = intronOrExonCanonical[0];
        data.totalExonIntronCanonical = intronOrExonCanonical[1];

        String[] intronOrExonWorst = getIntronOrExon(json, "worst");
        data.variantExonIntronWorst = intronOrExonWorst[0];
        data.totalExonIntronWorst = intronOrExonWorst[1];

        if (filters.spliceAiDsmax != null) {
            if (filters.spliceAiDsmax >= SpliceAIConnector.MAX_DS_UNLIKELY) {
                view.general.spliceAltering = Optional.ofNullable(getSpliceAltering(filters));
            }
        } else {
            view.general.spliceAltering = Optional.empty();
        }
    }

    private static String getSpliceAltering(AnfisaResultFilters filters) {
        return filters.spliceAltering;
    }

    private static String[] getIntronOrExon(JSONObject json, String kind) {
        List<String> introns = getFromTranscripts(json, "intron", kind);
        List<String> exons = getFromTranscripts(json, "exon", kind);
        List<String> e = (exons.size() > 0) ? exons : introns;
        if (e.size() == 0) {
            return new String[]{null, null};
        }

        List<String> index = new ArrayList<>();
        List<String> total = new ArrayList<>();
        for (String x : e) {
            if (x == null) continue;
            String[] split = x.split("/");
            index.add(split[0]);
            if (split.length > 1) {
                total.add(split[1]);
            }
        }

        return new String[]{
                String.join(",", index),
                String.join(",", total)
        };
    }


    private void createQualityTab(AnfisaResultView view, Variant variant, Map<String, Sample> samples) {
        if (variant == null || !(variant instanceof VariantVCF)) {
            return;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        CommonInfo commonInfo = variantContext.getCommonInfo();
        if (commonInfo == null) return;

        JSONObject q_all = new JSONObject();
        q_all.put("title", "All");
        q_all.put("strand_odds_ratio", toDouble(commonInfo.getAttribute("SOR")));
        q_all.put("mq", toDouble(commonInfo.getAttribute("MQ")));

        q_all.put("variant_call_quality", variantContext.getPhredScaledQual());
        q_all.put("qd", toDouble(commonInfo.getAttribute("QD")));
        q_all.put("fs", toDouble(commonInfo.getAttribute("FS")));
        q_all.put("ft", (commonInfo.isFiltered())
                ? commonInfo.getFilters().stream().collect(Collectors.toList())
                : Lists.newArrayList("PASS")
        );
        view.qualitySamples.add(q_all);

        String proband = SampleUtils.getProband(samples);
        if (proband == null) return;

        String mother = samples.get(proband).mother;
        String father = samples.get(proband).father;
        for (Map.Entry<String, Sample> entry : samples.entrySet()) {
            Sample sample = entry.getValue();
            String s = sample.id;
            JSONObject q_s = new JSONObject();
            if (s.equals(proband)) {
                q_s.put("title", String.format("Proband: %s", sample.name));
            } else if (s.equals(mother)) {
                q_s.put("title", String.format("Mother: %s", sample.name));
            } else if (s.equals(father)) {
                q_s.put("title", String.format("Father: %s", sample.name));
            } else {
                q_s.put("title", sample.name);
            }

            Genotype oGenotype = variantContext.getGenotype(sample.id);
            q_s.put("allelic_depth", oGenotype.getAnyAttribute("AD"));
            q_s.put("read_depth", oGenotype.getAnyAttribute("DP"));
            q_s.put("genotype_quality", getVariantGQ(variant, sample));
            view.qualitySamples.add(q_s);
        }
    }

    private void createGnomadTab(AnfisaExecuteContext context, String chromosome, Variant variant, Map<String, Sample> samples, JSONObject json, AnfisaResultFilters filters, AnfisaResultData data, AnfisaResultView view) {
        Double gnomadAf = context.gnomadAfFam;
        if (gnomadAf != null && Math.abs(gnomadAf) > 0.000001D) {
            for (String allele : alt_list(variant, samples, json)) {
                AnfisaResultView.GnomAD gnomAD = new AnfisaResultView.GnomAD();

                gnomAD.allele = allele;
                gnomAD.pli = getPLIByAllele(json, allele);
                gnomAD.proband = (isProbandHasAllele(variant, samples, allele)) ? "Yes" : "No";

                GnomadResult gnomadResult = getGnomadResult(variant, json, allele);
                if (gnomadResult != null) {
                    if (gnomadResult.exomes != null) {
                        gnomAD.exomeAn = gnomadResult.exomes.an;
                        gnomAD.exomeAf = gnomadResult.exomes.af;
                    }
                    if (gnomadResult.genomes != null) {
                        gnomAD.genomeAn = gnomadResult.genomes.an;
                        gnomAD.genomeAf = gnomadResult.genomes.af;
                    }

                    if (gnomadResult.overall != null) {
                        gnomAD.af = gnomadResult.overall.af;
                        gnomAD.hom = gnomadResult.overall.hom;
                        gnomAD.hem = gnomadResult.overall.hem;
                    }
                    if (gnomadResult.widePopmax != null) {
                        gnomAD.widePopmax = String.format(Locale.ENGLISH, "%s: %.5f [%s]",
                                gnomadResult.widePopmax.group.name(), gnomadResult.widePopmax.af, gnomadResult.widePopmax.an
                        );
                    }
                    gnomAD.url = gnomadResult.urls.stream().map(url -> url.toString()).toArray(String[]::new);
                }

                view.gnomAD.add(gnomAD);
            }
        } else {
            int p1 = lowest_coord(json) - 2;
            int p2 = highest_coord(json) + 1;

            AnfisaResultView.GnomAD gnomAD = new AnfisaResultView.GnomAD();
            gnomAD.url = new String[]{
                    String.format("https://gnomad.broadinstitute.org/region/%s-%s-%s", chromosome, p1, p2)
            };
            view.gnomAD.add(gnomAD);
        }
    }

    private static int lowest_coord(JSONObject json) {
        return Math.min(start(json), end(json));
    }

    private static int highest_coord(JSONObject json) {
        return Math.max(start(json), end(json));
    }

    private static int start(JSONObject json) {
        return json.getAsNumber("start").intValue();
    }

    private static int end(JSONObject json) {
        return json.getAsNumber("end").intValue();
    }

    private static List<Double> getPLIByAllele(JSONObject json, String allele) {
        List<JSONObject> transcripts = getTranscripts(json, "protein_coding");
        String key = "exacpli";
        List<Double> list = unique(
                transcripts.stream()
                        .filter(jsonObject -> jsonObject.containsKey(key))
                        .filter(jsonObject -> allele.equals(jsonObject.get("variant_allele")))
                        .map(jsonObject -> jsonObject.getAsNumber(key).doubleValue())
                        .collect(Collectors.toList())
        );
        if (list.size() > 0) {
            return list;
        }
        return null;
    }

    private static Set<String> getHGMDTags(Record record) {
        if (record.hgmdData == null)
            return new HashSet<>();
        return record.hgmdData.hgmdPmidRows.stream().map(hgmdPmidRow -> hgmdPmidRow.tag).collect(Collectors.toSet());
    }

    private void createDatabasesTab(JSONObject response, Record record, AnfisaResultData data, AnfisaResultView view) {
        if (data.hgmd != null) {
            view.databases.hgmd = data.hgmd;
            view.databases.hgmdHg38 = data.hgmdHg38;

            view.databases.hgmdTags = getHGMDTags(record).toArray(new String[0]);
            if (view.databases.hgmdTags.length == 0) view.databases.hgmdTags = null;

            view.databases.hgmdPhenotypes = record.hgmdData.phenotypes.toArray(new String[record.hgmdData.phenotypes.size()]);
            if (view.databases.hgmdPhenotypes.length == 0) view.databases.hgmdPhenotypes = null;

            view.databases.hgmdPmids = record.hgmdData.hgmdPmidRows.stream()
                    .map(hgmdPmidRow -> hgmdPmidRow.pmid).map(pmid -> linkToPmid(pmid)).toArray(String[]::new);
            if (view.databases.hgmdPmids.length == 0) view.databases.hgmdPmids = null;

            data.hgmdPmids = record.hgmdData.hgmdPmidRows.stream()
                    .map(hgmdPmidRow -> hgmdPmidRow.pmid).toArray(String[]::new);
            if (data.hgmdPmids.length == 0) data.hgmdPmids = null;
        } else {
            view.databases.hgmd = "Not Present";
        }
        view.databases.beaconUrl = data.beaconUrls;

        if (data.clinvarVariants != null) {
            view.databases.clinVarVariants = Arrays.stream(data.clinvarVariants).distinct().toArray(String[]::new);
        }
        if (data.clinvarSignificance != null) {
            view.databases.clinVarSignificance = Arrays.stream(data.clinvarSignificance).distinct().toArray(String[]::new);
        }
        if (data.clinvarPhenotypes != null) {
            view.databases.clinVarPhenotypes = Arrays.stream(data.clinvarPhenotypes).distinct().toArray(String[]::new);
        }
        for (String submitter : trustedSubmitters.keySet()) {
            view.databases.setField(String.format("%s_significance", submitter), data.getField(submitter));
        }
        view.databases.pubmedSearch = getTenwiseLink(response);
        view.databases.omim = getGenes(response).stream().map(gene ->
                String.format("https://omim.org/search/?search=approved_gene_symbol:%s&retrieve=geneMap", gene)
        ).toArray(String[]::new);
        view.databases.geneCards = getGenes(response).stream().map(gene ->
                String.format("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene)
        ).toArray(String[]::new);
    }

    private void createPredictionsTab(JSONObject json, AnfisaResultView view) {
        view.predictions.lofScore = getFromTranscripts(json, "loftool", "all")
                .stream().map(s -> Double.parseDouble(s)).sorted(Comparator.reverseOrder())
                .collect(Collectors.toList());
        view.predictions.lofScoreCanonical = getFromCanonicalTranscript(json, "loftool")
                .stream().map(s -> Double.parseDouble(s)).sorted(Comparator.reverseOrder())
                .collect(Collectors.toList());

        view.predictions.maxEntScan = getMaxEnt(json);

        view.predictions.polyphen = getFromTranscriptsList(json, "polyphen_prediction").stream().toArray(String[]::new);
        view.predictions.polyphen2Hvar = getFromTranscriptsList(json, "Polyphen2_HVAR_pred".toLowerCase()).stream().collect(Collectors.toList());
        view.predictions.polyphen2Hdiv = getFromTranscriptsList(json, "Polyphen2_HDIV_pred".toLowerCase()).stream().collect(Collectors.toList());
        view.predictions.polyphen2HvarScore = getFromTranscriptsList(json, "Polyphen2_HVAR_score".toLowerCase()).stream()
                .collect(Collectors.toList());
        view.predictions.polyphen2HdivScore = getFromTranscriptsList(json, "Polyphen2_HDIV_score".toLowerCase()).stream()
                .collect(Collectors.toList());
        view.predictions.sift = getFromTranscriptsList(json, "sift_prediction").stream().toArray(String[]::new);
        view.predictions.siftVEP = getFromTranscriptsList(json, "sift_pred").stream().toArray(String[]::new);
        view.predictions.siftScore = getFromTranscriptsList(json, "sift_score").stream().toArray(String[]::new);
        view.predictions.revel = getFromTranscriptsList(json, "revel_score").stream().map(s -> Double.parseDouble(s))
                .collect(Collectors.toList());
        view.predictions.mutationTaster = getFromTranscriptsList(json, "mutationtaster_pred").stream().toArray(String[]::new);
        view.predictions.fathmm = getFromTranscriptsList(json, "fathmm_pred").stream().toArray(String[]::new);
        view.predictions.caddPhred = getFromTranscriptsList(json, "cadd_phred").stream().map(s -> Double.parseDouble(s))
                .collect(Collectors.toList());
        view.predictions.caddRaw = getFromTranscriptsList(json, "cadd_raw").stream().map(s -> Double.parseDouble(s))
                .collect(Collectors.toList());
        view.predictions.mutationAssessor = getFromTranscriptsList(json, "mutationassessor_pred").stream().toArray(String[]::new);
    }

    private static List<String> getMaxEnt(JSONObject json) {
        List<JSONObject> transcripts = getTranscripts(json, "protein_coding");
        Set<String> x = new HashSet<>();
        for (JSONObject transcript : transcripts) {
            Number m1 = transcript.getAsNumber("maxentscan_ref");
            Number m2 = transcript.getAsNumber("maxentscan_alt");
            Number m3 = transcript.getAsNumber("maxentscan_diff");
            if (m1 != null && m2 != null && m3 != null) {
                String v = String.format("%s=%s-%s", m3, m1, m2);
                x.add(v);
            }
        }
        if (!x.isEmpty()) {
            return new ArrayList<String>(x);
        } else {
            return null;
        }
    }

    private void createBioinformaticsTab(GtfAnfisaResult gtfAnfisaResult, AnfisaExecuteContext anfisaExecuteContext, AnfisaResultData data, AnfisaResultView view) {
        AnfisaInput anfisaInput = anfisaExecuteContext.anfisaInput;
        Variant variant = anfisaExecuteContext.variant;
        JSONObject vepJson = anfisaExecuteContext.vepJson;

        view.bioinformatics.zygosity = getZygosity(vepJson, variant, anfisaInput.samples);
        view.bioinformatics.inheritedFrom = inherited_from(vepJson, variant, anfisaInput.samples);
        view.bioinformatics.distFromExonWorst = getDistanceFromExon(gtfAnfisaResult, vepJson, "worst");
        view.bioinformatics.distFromExonCanonical = getDistanceFromExon(gtfAnfisaResult, vepJson, "canonical");
        view.bioinformatics.conservation = buildConservation(anfisaExecuteContext);
        view.bioinformatics.speciesWithVariant = "";
        view.bioinformatics.speciesWithOthers = "";
        view.bioinformatics.maxEntScan = getMaxEnt(vepJson);
        view.bioinformatics.nnSplice = "";
        view.bioinformatics.humanSplicingFinder = "";
        view.bioinformatics.otherGenes = getOtherGenes(vepJson);
        view.bioinformatics.calledBy = getCallers(vepJson, variant, anfisaInput.samples).stream().toArray(String[]::new);
        view.bioinformatics.callerData = getCallersData(variant);
        view.bioinformatics.spliceAi = list_dsmax(data);

        String[] splice_ai_keys = new String[]{"AG", "AL", "DG", "DL"};
        if (!data.spliceAI.isEmpty()) {
            for (Map.Entry<String, SpliceAIResult.DictSql> entry : data.spliceAI.entrySet()) {
                for (String s : splice_ai_keys) {
                    String key = String.format("DS_%s", s);
                    float score = entry.getValue().getValue(key).floatValue();
                    if (score > 0.0f) {
                        String key2 = String.format("DP_%s", s);
                        int position = entry.getValue().getValue(key2).intValue();
                        String sPosition = String.valueOf(position);
                        if (position > 0) {
                            sPosition = "+" + sPosition;
                        }
                        view.bioinformatics.getSpliceAiValues(s).add(
                                String.format("%s: %s[%s]", entry.getKey(), String.format(Locale.ENGLISH, "%.4f", score), sPosition)
                        );
                    }
                }
            }
        }
    }

    public Conservation buildConservation(AnfisaExecuteContext context) {
        AnfisaInput anfisaInput = context.anfisaInput;
        Variant variant = context.variant;
        JSONObject vepJson = context.vepJson;
        Chromosome chromosome = variant.chromosome;
        String ref = getRef(variant, vepJson);
        List<String> alts = alt_list(variant, anfisaInput.samples, vepJson);
        if (alts.size() > 1) {
            return null;
        }
        String alt = alts.get(0);
        Position hgmdHG19 = getHg19Coordinates(context);
        Position hgmdHG38 = getHg38Coordinates(context);
        if (hgmdHG38 == null) {
            return null;
        }
        return conservationConnector.getConservation(chromosome, hgmdHG19, hgmdHG38, ref, alt);
    }

    private static Map<String, Float> list_dsmax(AnfisaResultData data) {
        Map<String, Float> result = new HashMap<>();
        if (data.spliceAI.isEmpty()) {
            return result;
        }
        result.put(
                "DS_AG",
                data.spliceAI.values().stream().map(dictSql -> dictSql.ds_ag).max(Float::compareTo).get()
        );
        result.put(
                "DS_AL",
                data.spliceAI.values().stream().map(dictSql -> dictSql.ds_al).max(Float::compareTo).get()
        );
        result.put(
                "DS_DG",
                data.spliceAI.values().stream().map(dictSql -> dictSql.ds_dg).max(Float::compareTo).get()
        );
        result.put(
                "DS_DL",
                data.spliceAI.values().stream().map(dictSql -> dictSql.ds_dl).max(Float::compareTo).get()
        );
        return result;
    }

    private static String[] getOtherGenes(JSONObject response) {
        Set<String> genes = new HashSet<>(getGenes(response));
        Set<String> allGenes = new HashSet<>(getFromTranscriptsByBiotype(response, "gene_symbol", "all"));

        Set<String> result = new HashSet<>();
        result.addAll(allGenes);
        result.removeAll(genes);
        return result.toArray(new String[result.size()]);
    }

    private String getStrHg38Coordinates(AnfisaExecuteContext context) {
        Chromosome chromosome = context.variant.chromosome;
        Position<Integer> positionHg38 = getHg38Coordinates(context);

        if (positionHg38 == null) {
            return String.format("%s:None", chromosome.toString());
        } else if (positionHg38.isSingle()) {
            return String.format("%s:%s", chromosome.toString(), positionHg38.start);
        } else {
            return String.format("%s:%s-%s",
                    chromosome.toString(), positionHg38.start, positionHg38.end
            );
        }
    }

    private Position<Integer> getHg38Coordinates(AnfisaExecuteContext context) {
        JSONObject vepJson = context.vepJson;
        Chromosome chromosome = context.variant.chromosome;

        return liftoverConnector.toHG38(chromosome, new Position<Long>(
                vepJson.getAsNumber("start").longValue(),
                vepJson.getAsNumber("end").longValue()
        ));
    }

    public ColorCode.Code getColorCode(JSONObject vepJson, AnfisaResultData data, Record record, AnfisaResultFilters filters) {
        String csq = data.mostSevereConsequence;
        Consequences.Severity msq = Consequences.severity(csq);

        ColorCode.Shape shape;
        if (msq == Consequences.Severity.DAMAGING)
            shape = ColorCode.Shape.CROSS;
        else
            shape = ColorCode.Shape.CIRCLE;

        Set<String> hgmdTags = getHGMDTags(record);
        hgmdTags.retainAll(ImmutableSet.of("DM", "DM?"));
        boolean hgmdDamaging = !hgmdTags.isEmpty();

        ColorCode.Color color;
        if (data.clinvarSignificance != null) {
            boolean benign = filters.clinvarTrustedBenign.orElse(filters.clinvarBenign);
            if (benign) {
                if (shape == ColorCode.Shape.CROSS || hgmdDamaging) {
                    color = ColorCode.Color.YELLOW;
                } else {
                    color = ColorCode.Color.GREEN;
                }
                return ColorCode.code(shape, color);
            }

            for (String s : data.clinvarSignificance) {
                s = s.toLowerCase();
                if (s.contains("pathogenic") && !s.contains("conflict")) {
                    color = ColorCode.Color.RED;
                } else if (s.contains("conflict") || s.contains("uncertain")) {
                    color = ColorCode.Color.YELLOW;
                } else {
                    continue;
                }

                return ColorCode.code(shape, color);
            }
        }
        if (hgmdDamaging) {
            return ColorCode.code(shape, ColorCode.Color.RED);
        }

        int best = 100;
        int worst = 0;
        for (String tool : ColorCode.allInSilicoTools()) {
            List<String> rawValues = getFromTranscriptsList(vepJson, tool);
            for (String rawValue : rawValues) {
                int value = ColorCode.inSilicoPrediction(tool, rawValue);
                if (value == 0)
                    continue;
                if (value > worst)
                    worst = value;
                if (value < best)
                    best = value;
            }
        }

        if (worst >= 30)
            color = ColorCode.Color.RED;
        else if (best <= 10)
            color = ColorCode.Color.GREEN;
        else if (best < 100)
            color = ColorCode.Color.YELLOW;
        else if (shape == ColorCode.Shape.CROSS)
            return ColorCode.code(shape, ColorCode.Color.YELLOW);
        else
            return null;

        return ColorCode.code(shape, color);
    }

    private String[] getTenwiseLink(JSONObject response) {
        List<String> hgncIds = getHgncIds(response);
        return hgncIds.stream().map(hgncId ->
                String.format("https://www.tenwiseapps.nl/publicdl/variant_report/HGNC_%s_variant_report.html", hgncId)
        ).toArray(String[]::new);
    }

    public String altString(Variant variant, Map<String, Sample> samples, JSONObject response) {
        return String.join(",", alt_list(variant, samples, response));
    }

    private static List<String> alt_list(Variant variant, Map<String, Sample> samples, JSONObject vepJson) {
        if (variant != null && variant instanceof VariantVCF) {
            VariantContext variantContext = ((VariantVCF) variant).variantContext;
            List<String> alleles = variantContext.getAlleles()
                    .stream().map(allele -> allele.getBaseString()).collect(Collectors.toList());
            List<String> alt_allels = variantContext.getAlternateAlleles()
                    .stream().map(allele -> allele.getBaseString()).collect(Collectors.toList());
            if (samples != null) {
                Map<String, Long> counts = new HashMap<>();
                for (Genotype genotype : variantContext.getGenotypes()) {
                    int[] ad = genotype.getAD();
                    if (ad == null || ad.length == 0) {
                        return alt_allels;
                    }
                    for (int i = 0; i < alleles.size(); i++) {
                        String al = alleles.get(i);
                        long n = ad[i];
                        counts.put(al, counts.getOrDefault(al, 0L) + n);
                    }
                }
                List<String> tmp_alt_allels = alt_allels.stream()
                        .filter(s -> counts.containsKey(s) && counts.get(s) > 0)
                        .collect(Collectors.toList());
                if (!tmp_alt_allels.isEmpty()) {
                    return tmp_alt_allels;
                } else {
                    return alt_allels;
                }
            } else {
                return alt_allels;
            }
        } else {
            return getAlts1(vepJson);
        }
    }

    private static List<String> getAlts1(JSONObject vepJson) {
        String[] ss = getAllele(vepJson).split("/");
        List<String> result = new ArrayList<>();
        for (int i = 1; i < ss.length; i++) {
            result.add(ss[i]);
        }
        return result;
    }

    private static String getAllele(JSONObject vepJson) {
        return vepJson.getAsString("allele_string");
    }

    private static String getRef(Variant variant, JSONObject vepJson) {
        if (variant != null && variant instanceof VariantVCF) {
            VariantContext variantContext = ((VariantVCF) variant).variantContext;
            return variantContext.getReference().getBaseString();
        } else {
            return getRef1(vepJson);
        }
    }

    private static String getRef1(JSONObject response) {
        return response.getAsString("allele_string").split("/")[0];
    }

    private static boolean isProbandHasAllele(Variant variant, Map<String, Sample> samples, String alt) {
        if (samples == null || variant == null) {
            return false;
        }

        String probandGenotype = (String) getGenotypes(variant, samples)[0];
        if (probandGenotype == null) {
            return false;
        }
        Set<String> set1 = Arrays.stream(probandGenotype.split("/")).collect(Collectors.toSet());
        return set1.contains(alt);
    }

    private static Integer getMinGQ(Variant variant, Map<String, Sample> samples) {
        if (variant == null || samples == null) {
            return null;
        }

        Integer GQ = null;
        for (Sample s : samples.values()) {
            Integer gq = getVariantGQ(variant, s);
            if (gq != null && gq != 0) {
                if (GQ == null || gq < GQ) {
                    GQ = gq;
                }
            }
        }
        return GQ;
    }

    private static Integer getProbandGQ(Variant variant, Map<String, Sample> samples) {
        if (samples == null || variant == null) {
            return null;
        }
        return getVariantGQ(variant, samples.get(SampleUtils.getProband(samples)));
    }

    private static Integer getVariantGQ(Variant variant, Sample s) {
        if (variant == null || !(variant instanceof VariantVCF)) {
            return null;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        int valie = variantContext.getGenotype(s.id).getGQ();
        return (valie != -1) ? valie : null;
    }

    private Long getSeverity(JSONObject response) {
        String csq = response.getAsString("most_severe_consequence");
        int n = AnfisaVariant.SEVERITY.size();
        for (int s = 0; s < n; s++) {
            if (AnfisaVariant.SEVERITY.get(s).contains(csq)) {
                return Long.valueOf(n - s - 2);
            }
        }
        return null;
    }

    private static List<Object> getDistanceFromExon(GtfAnfisaResult gtfAnfisaResult, JSONObject json, String kind) {
        List<Object[]> distFromBoundary;
        if ("canonical".equals(kind)) {
            distFromBoundary = gtfAnfisaResult.distFromBoundaryCanonical;
        } else if ("worst".equals(kind)) {
            distFromBoundary = gtfAnfisaResult.distFromBoundaryWorst;
        } else {
            throw new RuntimeException("Unknown kind: " + kind);
        }
        if (distFromBoundary != null) {
            return unique(
                    distFromBoundary.stream().map(objects -> (Long) objects[0]).collect(Collectors.toList())
            );
        }

        return unique(
                getHgvsList(json, "c", kind).stream()
                        .map(hgvcs -> getDistanceHgvsc(hgvcs))
                        .collect(Collectors.toList()),
                "Exonic"
        );
    }

    private static List<Character> hgvs_signs = Lists.newArrayList('-', '+', '*');

    private static Long getDistanceHgvsc(String hgvcs) {
        String[] sChunks = hgvcs.split(":");
        String coord = null;
        for (int i = 1; i < sChunks.length; i++) {
            String chunk = sChunks[i];
            char ch = chunk.charAt(0);
            if (ch == 'c' || ch == 'p') {
                coord = chunk;
                break;
            }
        }
        if (coord == null) {
            return null;
        }
        String[] xx = coord.split("\\.");
        Integer d = null;
        try {
            for (String x : xx[1].split("_")) {
                Integer sign = null;
                Integer p1 = null;
                Integer p2 = null;
                while (!Character.isDigit(x.charAt(0)) && !hgvs_signs.contains(x.charAt(0))) {
                    x = x.substring(1);
                }
                int end = x.length();
                for (int i = 0; i < end; i++) {
                    char c = x.charAt(i);
                    if (Character.isDigit(c)) {
                        continue;
                    }
                    if (hgvs_signs.contains(c)) {
                        int p0 = (sign == null) ? 0 : sign + 1;
                        p1 = (i > p0) ? Integer.parseInt(x.substring(p0, i)) : 0;
                        sign = i;
                    }
                    if (Character.isAlphabetic(c)) {
                        end = i;
                        break;
                    }
                }
                if (p1 != null && sign != null) {
                    p2 = Integer.parseInt(x.substring(sign + 1, end));
                }
                if (p2 != null) {
                    if (d == null || d > p2) {
                        d = p2;
                    }
                }
            }
            return (d != null) ? d.longValue() : null;
        } catch (Exception e) {
            log.error("Exception ", e);
            return null;
        }
    }

    public String getLabel(AnfisaExecuteContext context) {
        List<String> genes = getGenes(context.vepJson);
        String gene;
        if (genes.size() == 0) {
            gene = "None";
        } else if (genes.size() < 3) {
            gene = String.join(",", genes);
        } else {
            gene = "...";
        }

        String vstr = str(context);

        return String.format("[%s] %s", gene, vstr);
    }

    private static List<String> getGenes(JSONObject response) {
        return getFromTranscriptsList(response, "gene_symbol");
    }

    public List<String> getHgncIds(JSONObject response) {
        return getFromTranscriptsList(response, "hgnc_id");
    }

    /**
     * def get_most_severe_transcripts(self):
     * msq = self.get_msq()
     * return [t for t in self.get_transcripts() if (msq in t.get("consequence_terms"))]
     */
    private static List<JSONObject> getMostSevereTranscripts(JSONObject json) {
        String msq = getMsq(json);
        return getTranscripts(json, "protein_coding").stream()
                .filter(jsonObject -> {
                    JSONArray consequenceTerms = (JSONArray) jsonObject.get("consequence_terms");
                    return (consequenceTerms != null && consequenceTerms.contains(msq));
                })
                .collect(Collectors.toList());
    }

    private static List<JSONObject> getCanonicalTranscripts(JSONObject json) {
        return getTranscripts(json, "protein_coding").stream()
                .filter(jsonObject -> jsonObject.containsKey("canonical"))
                .collect(Collectors.toList());
    }


    private static List<String> getFromTranscriptsList(JSONObject json, String key) {
        return getFromTranscriptsByBiotype(json, key, "protein_coding");
    }

    private static List<JSONObject> getTranscripts(JSONObject response, String biotype) {
        List<JSONObject> result = new ArrayList<>();
        JSONArray jTranscriptConsequences = (JSONArray) response.get("transcript_consequences");
        if (jTranscriptConsequences != null) {
            for (Object oItem : jTranscriptConsequences) {
                JSONObject item = (JSONObject) oItem;
                if (biotype == null || biotype.toUpperCase().equals("ALL")) {
                    result.add(item);
                } else if (item.get("biotype").equals(biotype)) {
                    result.add(item);
                }
            }
        }
        return result;
    }

    private static List<String> getFromTranscriptsByBiotype(JSONObject response, String key, String biotype) {
        List<String> result = new ArrayList<>();

        for (JSONObject item : getTranscripts(response, biotype)) {
            Object oValue = item.get(key);
            if (oValue == null) continue;
            if (oValue instanceof JSONArray) {
                for (Object io : (JSONArray) oValue) {
                    String i = String.valueOf(io);
                    if (!result.contains(i)) {
                        result.add(i);
                    } else {
                        //TODO Ulitin V. Необходимо выяснить, нужна ли такая особенность пострения уникального списка
                        //Дело в том, что в python реализации косвенно получался такой результат,
                        //непонятно, это сделано специально или нет, если не важна сортировака, то заменить на обычный Set
                        result.remove(i);
                        result.add(0, i);
                    }
                }
            } else {
                String value = String.valueOf(oValue);
                if (!result.contains(value)) {
                    result.add(value);
                } else {
                    //TODO Ulitin V. Необходимо выяснить, нужна ли такая особенность пострения уникального списка
                    //Дело в том, что в python реализации косвенно получался такой результат,
                    //непонятно, это сделано специально или нет, если не важна сортировака, то заменить на обычный Set
                    result.remove(value);
                    result.add(0, value);
                }
            }
        }
        return result;
    }

    private static List<String> getFromWorstTranscript(JSONObject json, String key) {
        return unique(
                getMostSevereTranscripts(json).stream()
                        .filter(jsonObject -> jsonObject.containsKey(key))
                        .map(jsonObject -> jsonObject.getAsString(key))
                        .collect(Collectors.toList())
        );
    }

    private static List<String> getFromCanonicalTranscript(JSONObject json, String key) {
        return unique(
                getCanonicalTranscripts(json).stream()
                        .filter(jsonObject -> jsonObject.containsKey(key))
                        .flatMap(jsonObject -> {
                            Object value = jsonObject.get(key);
                            if (value instanceof String) {
                                return Stream.of((String) value);
                            } else if (value instanceof Number) {
                                return Stream.of(String.valueOf((Number) value));
                            } else {
                                return ((JSONArray) value)
                                        .stream()
                                        .map(o -> (String) o);
                            }
                        })
                        .collect(Collectors.toList())
        );
    }

    private static List<String> getFromTranscripts(List<JSONObject> transcripts, String key, String source) {
        return unique(
                transcripts.stream()
                        .filter(jsonObject -> source.equals(jsonObject.getAsString("source")))
                        .map(jsonObject -> jsonObject.getAsString(key))
                        .collect(Collectors.toList())
        );
    }

    private static List<String> getFromTranscripts(JSONObject json, String key, String type) {
        if ("all".equals(type)) {
            return getFromTranscriptsList(json, key);
        } else if ("canonical".equals(type)) {
            return getFromCanonicalTranscript(json, key);
        } else if ("worst".equals(type)) {
            return getFromWorstTranscript(json, key);
        } else {
            throw new RuntimeException("Unknown type: " + type);
        }
    }

    private static List<String> getHgvsList(JSONObject json, String type, String kind) {
        if ("c".equals(type)) {
            return getFromTranscripts(json, "hgvsc", kind);
        } else if ("p".equals(type)) {
            return getFromTranscripts(json, "hgvsp", kind);
        } else {
            List<String> result = new ArrayList<>();
            result.addAll(getFromTranscripts(json, "hgvsc", kind));
            result.addAll(getFromTranscripts(json, "hgvsp", kind));
            return result;
        }
    }

    private static String hgvcsPos(String str, String type, boolean withPattern) {
        String pattern = String.format(":%s.", type);
        if (!str.contains(pattern)) {
            return null;
        }
        String x = str.split(pattern)[1];
        if ("p".equals(type)) {
            x = convert_p(x);
        }

        if (withPattern) {
            return String.format("%s%s", pattern.substring(1), x);
        } else {
            return x;
        }
    }

    private static String convert_p(String x) {
        List<Character> protein1 = new ArrayList<>();
        List<Character> pos = new ArrayList<>();
        List<Character> protein2 = new ArrayList<>();
        int state = 0;
        for (char c : x.toCharArray()) {
            if (state == 0) {
                if (Character.isLetter(c)) {
                    protein1.add(c);
                    continue;
                }
                state = 2;
            }
            if (state == 2) {
                if (Character.isDigit(c)) {
                    pos.add(c);
                    continue;
                }
                state = 3;
            }
            if (state == 3) {
                protein2.add(c);
            } else {
                break;
            }
        }
        String p1 = protein1.stream().map(c -> c.toString()).collect(Collectors.joining());
        String p2 = protein2.stream().map(c -> c.toString()).collect(Collectors.joining());
        String rpos = pos.stream().map(c -> c.toString()).collect(Collectors.joining());
        String rprotein1 = proteins_3_to_1.getOrDefault(p1, p1);
        String rprotein2 = proteins_3_to_1.getOrDefault(p2, p2);
        return String.format("%s%s%s", rprotein1, rpos, rprotein2);
    }

    private static List<String> getPos(JSONObject json, String type, String kind) {
        List<String> hgvsList = getHgvsList(json, type, kind);
        List<String> poss = hgvsList.stream()
                .map(hgvcs -> hgvcsPos(hgvcs, type, true))
                .collect(Collectors.toList());
        return unique(poss);
    }

    private static List<String>[] getPosTpl(JSONObject json, String type) {
        Set<String> ss = new HashSet<>();

        List<String> c_worst = getPos(json, type, "worst");
        ss.addAll(c_worst);

        List<String> c_canonical = getPos(json, type, "canonical");
        ss.addAll(c_canonical);

        List<String> c_other = getPos(json, type, "all");
        if (c_other.isEmpty()) {
            c_other = Collections.emptyList();
        } else {
            ss.removeAll(c_other);
            c_other = (ss.isEmpty()) ? null : new ArrayList<>(ss);
        }

        return new List[]{
                c_worst, c_canonical, c_other
        };
    }

    public String str(AnfisaExecuteContext context) {
        AnfisaInput anfisaInput = context.anfisaInput;
        String str = getStrHg19Coordinates(context);
        if (isSnv(context.vepJson)) {
            return String.format("%s  %s>%s",
                    str,
                    ref(context.vepJson, context.variant),
                    altString(context.variant, anfisaInput.samples, context.vepJson)
            );
        } else {
            String variantClass = getVariantClass(context.vepJson);
            return String.format("%s %s", str, (variantClass != null) ? variantClass : "None");
        }
    }

    public Position getHg19Coordinates(AnfisaExecuteContext context) {
        JSONObject vepJson = context.vepJson;
        return new Position(
                vepJson.getAsNumber("start").longValue(),
                vepJson.getAsNumber("end").longValue()
        );
    }

    public String getStrHg19Coordinates(AnfisaExecuteContext context) {
        Position hg19Coordinates = getHg19Coordinates(context);
        return vstr(
                context.variant.chromosome.getChromosome(),
                hg19Coordinates
        );
    }
//
//    public String vstr(String c, long s, long e) {
//        if (s == e) {
//            return String.format("%s:%s", c, s);
//        } else {
//            return String.format("%s:%s-%s", c, s, e);
//        }
//    }

    public String vstr(String c, Position hg19Coordinates) {
        if (hg19Coordinates.isSingle()) {
            return String.format("%s:%s", c, hg19Coordinates.start);
        } else {
            return String.format("%s:%s-%s", c, hg19Coordinates.start, hg19Coordinates.end);
        }
    }

    private static String getMsq(JSONObject json) {
        return json.getAsString("most_severe_consequence");
    }

    public boolean isSnv(JSONObject json) {
        return "SNV".equals(getVariantClass(json));
    }

    private String getVariantClass(JSONObject json) {
        return json.getAsString("variant_class");
    }

    public String linkToPmid(String pmid) {
        return String.format("https://www.ncbi.nlm.nih.gov/pubmed/%s", pmid);
    }

    private static String encodeToAscii(String s) {
        String regex = "[\\p{InCombiningDiacriticalMarks}\\p{IsLm}\\p{IsSk}]+";
        try {
            return new String(s.replaceAll(regex, "").getBytes("ascii"), "ascii");
        } catch (UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
    }

    private static String ref(JSONObject json, Variant variant) {
        if (variant != null && variant instanceof VariantVCF) {
            VariantContext variantContext = ((VariantVCF) variant).variantContext;
            return variantContext.getReference().getBaseString();
        } else {
            return ref1(json);
        }
    }

    private static String ref1(JSONObject json) {
        String s = json.getAsString("allele_string");
        return s.split("/")[0];
    }

    private static Object[] getGenotypes(Variant variant, Map<String, Sample> samples) {
        String empty = "Can not be determined";
        String proband = SampleUtils.getProband(samples);
        if (proband == null) {
            return new Object[]{null, null, null, null};
        }
        String probandGenotype = getGtBasesGenotype(variant, proband);
        if (probandGenotype == null) {
            probandGenotype = empty;
        }

        String mother = samples.get(proband).mother;
        if ("0".equals(mother)) {
            mother = null;
        }
        String maternalGenotype = (mother != null) ? getGtBasesGenotype(variant, mother) : null;
        if (mother != null && maternalGenotype == null) {
            maternalGenotype = empty;
        }

        String father = samples.get(proband).father;
        if ("0".equals(father)) {
            father = null;
        }
        String paternalGenotype = (father != null) ? getGtBasesGenotype(variant, father) : null;
        if (father != null && paternalGenotype == null) {
            paternalGenotype = empty;
        }

        String finalProbandGenotype = probandGenotype;
        String finalMaternalGenotype = maternalGenotype;
        String finalPaternalGenotype = paternalGenotype;
        List<String> otherGenotypes = samples.keySet().stream()
                .map(genotype -> getGtBasesGenotype(variant, genotype))
                .filter(gtBases -> gtBases != null)
                .filter(gtBases -> !gtBases.equals(finalProbandGenotype))
                .filter(gtBases -> !gtBases.equals(finalMaternalGenotype))
                .filter(gtBases -> !gtBases.equals(finalPaternalGenotype))
                .distinct()
                .collect(Collectors.toList());

        return new Object[]{probandGenotype, maternalGenotype, paternalGenotype, otherGenotypes};
    }

    private static String getGtBasesGenotype(Variant variant, String genotype) {
        if (variant == null || !(variant instanceof VariantVCF)) {
            return null;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        Genotype oGenotype = variantContext.getGenotype(genotype);
        if (oGenotype.isCalled()) {
            return oGenotype.getGenotypeString();
        } else {
            return null;
        }
    }

    private static String getMostSevere(List<String> consequenceTerms) {
        for (String item : AnfisaVariant.CONSEQUENCES) {
            if (consequenceTerms.contains(item)) {
                return item;
            }
        }
        return null;
    }

    private static LinkedHashSet<String> getRawCallers(VariantContext variantContext) {
        LinkedHashSet<String> callers = new LinkedHashSet<>();
        CommonInfo commonInfo = variantContext.getCommonInfo();
        for (String caller : AnfisaVariant.CALLERS) {
            if (commonInfo.hasAttribute(caller)) {
                if (AnfisaVariant.BGM_BAYES_DE_NOVO.equals(caller) &&
                        Double.parseDouble(commonInfo.getAttribute(caller).toString()) < 0
                ) {
                    //Отрицательное число, означает, что при работе коллера произошла ошибка…
                    callers.add(AnfisaVariant.BGM_BAYES_DE_NOVO_S1);
                    continue;
                }
                callers.add(caller);
            }
        }
        return callers;
    }


    private static LinkedHashSet<String> getCallers(JSONObject json, Variant variant, Map<String, Sample> samples) {
        if (samples == null || variant == null || !(variant instanceof VariantVCF)) {
            return new LinkedHashSet();
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        LinkedHashSet<String> callers = getRawCallers(variantContext);

        Object[] genotypes = getGenotypes(variant, samples);
        String probandGenotype = (String) genotypes[0];
        String maternalGenotype = (String) genotypes[1];
        String paternalGenotype = (String) genotypes[2];
        if (probandGenotype == null || maternalGenotype == null || paternalGenotype == null) {
            return callers;
        }

        String ref = ref(json, variant);
        List<String> alt_set = alt_list(variant, samples, json);
        Set<String> p_set = Arrays.stream(probandGenotype.split("/")).collect(Collectors.toSet());
        Set<String> m_set = Arrays.stream(maternalGenotype.split("/")).collect(Collectors.toSet());
        Set<String> f_set = Arrays.stream(paternalGenotype.split("/")).collect(Collectors.toSet());

        for (String alt : alt_set) {
            if (p_set.contains(alt) && !m_set.contains(alt) && !f_set.contains(alt)) {
                callers.add("GATK_DE_NOVO");
                break;
            }
        }

        if (p_set.size() == 1 && !Collections.disjoint(alt_set, p_set)) {
            if (m_set.size() == 2 && f_set.size() == 2 && !Collections.disjoint(alt_set, m_set) && !Collections.disjoint(alt_set, f_set)) {
                callers.add("GATK_HOMO_REC");
            }
        }

        if (p_set.size() == 1 && p_set.contains(ref)) {
            if (m_set.size() == 2 && f_set.size() == 2 && m_set.contains(ref) && f_set.contains(ref)) {
                callers.add("GATK_HOMOZYGOUS");
            }
        }

        if (callers.isEmpty()) {
            String inheritance = inherited_from(json, variant, samples);
            if ("De-Novo".equals(inheritance)) {
                throw new RuntimeException("Inconsistent inheritance");
            }
            if (!"Inconclusive".equals(inheritance)) {
                callers.add(String.format("INHERITED_FROM: %s", inheritance));
            }
        }

        return callers;
    }

    private static Map<String, Serializable> getCallersData(Variant variant) {
        if (variant == null || !(variant instanceof VariantVCF)) {
            return Collections.emptyMap();
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        CommonInfo commonInfo = variantContext.getCommonInfo();

        Map<String, Serializable> result = new HashMap<>();
        LinkedHashSet<String> callers = getRawCallers(variantContext);
        for (String c : callers) {
            if (commonInfo.hasAttribute(c)) {
                Object value = commonInfo.getAttribute(c);
                if (value instanceof Boolean) {
                    result.put(c, (Boolean) value);
                } else if (value instanceof String) {
                    //В первой реализации было приведение к Long, но как оказалось,
                    //в файле bgm9001_wgs_xbrowse.vep.vep.json есть значения вещественного типа
                    result.put(c, Double.parseDouble((String) value));
                } else {
                    throw new RuntimeException("Not support!");
                }
            }
        }
        return result;
    }

    private static String getZygosity(JSONObject json, Variant variant, Map<String, Sample> samples) {
        if (samples == null || variant == null || !(variant instanceof VariantVCF)) {
            return null;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        String chr = variant.chromosome.getChar();

        String genotype = (String) getGenotypes(variant, samples)[0];
        if (genotype == null) {
            return null;
        }
        List<String> set1 = Arrays.stream(genotype.split("/")).distinct().collect(Collectors.toList());

        if ("X".equals(chr.toUpperCase()) && proband_sex(variantContext, samples) == 1) {
            return "X-linked";
        }
        if (set1.size() == 1) {
            return "Homozygous";
        }
        if (set1.size() == 2) {
            return "Heterozygous";
        }
        return "Unknown";
    }

    private static String inherited_from(JSONObject json, Variant variant, Map<String, Sample> samples) {
        if (samples == null || variant == null || !(variant instanceof VariantVCF)) {
            return null;
        }
        VariantContext variantContext = ((VariantVCF) variant).variantContext;

        Object[] genotypes = getGenotypes(variant, samples);
        String probandGenotype = (String) genotypes[0];
        String maternalGenotype = (String) genotypes[1];
        String paternalGenotype = (String) genotypes[2];

        String chr = variant.chromosome.getChar();

        if ("X".equals(chr.toUpperCase()) && proband_sex(variantContext, samples) == 1) {
            if (probandGenotype.equals(maternalGenotype)) {
                return "Mother";
            } else {
                return "Inconclusive";
            }
        }

        if (maternalGenotype != null && maternalGenotype.equals(paternalGenotype)) {
            if (probandGenotype.equals(maternalGenotype)) {
                return "Both parents";
            } else {
                return "Both parents";
            }
        }

        if (probandGenotype != null && probandGenotype.equals(maternalGenotype)) {
            return "Mother";
        }
        if (probandGenotype != null && probandGenotype.equals(paternalGenotype)) {
            return "Father";
        }
        return "Inconclusive";
    }

    private static Integer proband_sex(VariantContext variantContext, Map<String, Sample> samples) {
        if (samples == null || variantContext == null) {
            return null;
        }
        String proband = SampleUtils.getProband(samples);
        return samples.get(proband).sex;
    }

    private static <T> List<T> unique(List<T> lst) {
        return unique(lst, null);
    }

    /**
     * TODO Ulitin V. Необходимо выяснить, нужна ли такая особенность пострения уникального списка
     * Дело в том, что в python реализации метод основан на set - в итоге проблема с сохранением порядка,
     * непонятно, это сделано специально или нет, если не важна сортировака,
     * то необходимо избавится от этого метода и перейти на Set
     *
     * @param lst
     * @return
     */
    private static <T> List<T> unique(List<T> lst, T replace_None) {
        if (lst == null) {
            return lst;
        }
        List<T> result = new ArrayList<>();
        for (T item : lst) {
            if (!result.contains(item)) {
                result.add(item);
            } else {
                result.remove(item);
                result.add(0, item);
            }
        }
        if (replace_None != null) {
            if (result.contains(null)) {
                result.remove(null);
                if (!result.contains(replace_None)) {
                    result.add(replace_None);
                }
            }
        }
        return result;
    }

    private static Double toDouble(Object value) {
        if (value == null) return null;
        if (value instanceof Number) {
            return ((Number) value).doubleValue();
        } else if (value instanceof String) {
            return Double.parseDouble((String) value);
        } else {
            throw new RuntimeException("Not support type");
        }
    }

    private static double toPrimitiveDouble(Object value) {
        Double p = toDouble(value);
        if (p == null) {
            return 0;
        } else {
            return p;
        }
    }

    @Override
    public void close() {

    }
}
