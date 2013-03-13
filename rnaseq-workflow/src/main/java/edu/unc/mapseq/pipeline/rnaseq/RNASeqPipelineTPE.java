package edu.unc.mapseq.pipeline.rnaseq;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class RNASeqPipelineTPE extends ThreadPoolExecutor {

    public RNASeqPipelineTPE() {
        super(20, 20, 5L, TimeUnit.MINUTES, new LinkedBlockingQueue<Runnable>());
    }

}