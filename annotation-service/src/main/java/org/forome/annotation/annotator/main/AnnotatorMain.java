package org.forome.annotation.annotator.main;

import org.forome.annotation.annotator.main.argument.*;
import org.forome.annotation.inventory.Inventory;
import org.forome.annotation.utils.AppVersion;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * cd /data/bgm/cases/bgm9001/
 * java -cp /home/vulitin/deploy/annotationservice/exec/annotation.jar org.forome.annotation.annotator.main.AnnotatorMain -config /home/vulitin/deploy/annotationservice/exec/config.json -vcf bgm9001_wgs_xbrowse.vep.vcf -vepjson bgm9001_wgs_xbrowse.vep.vep.json -output bgm9001_wgs_xbrowse.out.json
 * Для 6 милионов 37:09:11.460
 */
public class AnnotatorMain {

    private final static Logger log = LoggerFactory.getLogger(AnnotatorMain.class);

    public static void main(String[] args) {
        Arguments arguments;
        try {
            ParserArgument argumentParser = new ParserArgument(args);
            arguments = argumentParser.arguments;
        } catch (Throwable e) {
            log.error("Exception arguments parser", e);
            System.exit(2);
            return;
        }

        if (arguments instanceof ArgumentsVersion) {
            System.out.println("Version: " + AppVersion.getVersion());
            System.out.println("Version Format: " + AppVersion.getVersionFormat());
        } else if (arguments instanceof ArgumentsInventory) {
            ArgumentsInventory argumentsInventory = (ArgumentsInventory) arguments;
            Inventory inventory = new Inventory.Builder(argumentsInventory.pathInventory).build();
            AnnotationConsole annotationConsole = new AnnotationConsole(
                    argumentsInventory.config,
                    inventory.caseName, inventory.casePlatform,
                    inventory.famFile, inventory.patientIdsFile,
                    inventory.vcfFile, inventory.vepJsonFile,
                    0,
                    inventory.outFile
            );
            annotationConsole.execute();
        } else if (arguments instanceof ArgumentsAnnotation) {
            ArgumentsAnnotation argumentsAnnotation = (ArgumentsAnnotation) arguments;
            AnnotationConsole annotationConsole = new AnnotationConsole(
                    argumentsAnnotation.config,
                    argumentsAnnotation.caseName, argumentsAnnotation.casePlatform,
                    argumentsAnnotation.pathFam, argumentsAnnotation.patientIdsFile,
                    argumentsAnnotation.pathVcf, argumentsAnnotation.pathVepJson,
                    argumentsAnnotation.start,
                    argumentsAnnotation.pathOutput
            );
            annotationConsole.execute();
        } else {
            log.error("Unknown arguments");
            System.exit(3);
            return;
        }
    }

}
