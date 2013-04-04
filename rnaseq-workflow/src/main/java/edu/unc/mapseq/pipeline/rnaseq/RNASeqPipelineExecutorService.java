package edu.unc.mapseq.pipeline.rnaseq;

import java.util.Timer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class RNASeqPipelineExecutorService {

    private final Logger logger = LoggerFactory.getLogger(RNASeqPipelineExecutorService.class);

    private final Timer mainTimer = new Timer();

    private RNASeqPipelineBeanService pipelineBeanService;

    public void start() throws Exception {
        logger.info("ENTERING stop()");

        long delay = 1 * 60 * 1000;
        long period = 5 * 60 * 1000;

        RNASeqPipelineExecutorTask task = new RNASeqPipelineExecutorTask();
        task.setPipelineBeanService(pipelineBeanService);
        mainTimer.scheduleAtFixedRate(task, delay, period);

    }

    public void stop() throws Exception {
        logger.info("ENTERING stop()");
        mainTimer.purge();
        mainTimer.cancel();
    }

    public RNASeqPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(RNASeqPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

}
