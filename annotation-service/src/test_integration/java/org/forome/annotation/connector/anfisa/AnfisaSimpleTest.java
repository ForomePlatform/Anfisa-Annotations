package org.forome.annotation.connector.anfisa;

import org.forome.annotation.AnfisaBaseTest;
import org.forome.annotation.annotator.AnnotatorTest;
import org.forome.annotation.connector.anfisa.struct.AnfisaResult;
import org.forome.annotation.controller.GetAnfisaJSONController;
import org.forome.annotation.struct.Chromosome;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.concurrent.ExecutionException;

public class AnfisaSimpleTest extends AnfisaBaseTest {

    private final static Logger log = LoggerFactory.getLogger(AnnotatorTest.class);

    @Test
    public void test() throws Exception {
        GetAnfisaJSONController.build(
                request(new Chromosome("1"), 6484880, 6484880, "G")
        ).toJSONString();
//		anfisaConnector.request("1", 6500660, 6500660,"A").get();
//		anfisaConnector.request("1", 6501044, 6501044,"G").get();
//      anfisaConnector.request("1", 12040324, 12040324,"G").get();
//		anfisaConnector.request("1", 12065841, 12065841,"T").get();
//		anfisaConnector.request("1", 16351275, 16351275,"G").get();

//		anfisaConnector.request("1", 12065841, 12065841,"T").get().forEach(anfisaResult -> {
//			GetAnfisaJSONController.build(anfisaResult).toJSONString();
//		});
//        anfisaConnector.request("1", 6484880, 6484880, "G").get().forEach(anfisaResult -> {
//            GetAnfisaJSONController.build(anfisaResult).toJSONString();
//        });
//        anfisaConnector.request("1", 33475993, 33475993, "GC").get().forEach(anfisaResult -> {
//            GetAnfisaJSONController.build(anfisaResult).toJSONString();
//        });
//        anfisaConnector.request("1", 12058802, 12058802, "C").get().forEach(anfisaResult -> {
//            GetAnfisaJSONController.build(anfisaResult).toJSONString();
//        });
//        anfisaConnector.request("1", 16351275, 16351275, "G").get().forEach(anfisaResult -> {
//            GetAnfisaJSONController.build(anfisaResult).toJSONString();
//        });
//        anfisaConnector.request("1", 12040324, 12040324, "G").get().forEach(anfisaResult -> {
//            GetAnfisaJSONController.build(anfisaResult).toJSONString();
//        });


//        log.info("{}", GetAnfisaJSONController.build(
//                request(new Chromosome("1"), 216062306, 216062306, "G")
//                ).toJSONString()
//        );
    }

    private AnfisaResult request(Chromosome chromosome, long start, long end, String alternative) throws ExecutionException, InterruptedException {
        AnfisaResult anfisaResult = anfisaConnector.request(new Chromosome("1"), 216062306, 216062306, "G").get();
        return anfisaResult;
    }
}
