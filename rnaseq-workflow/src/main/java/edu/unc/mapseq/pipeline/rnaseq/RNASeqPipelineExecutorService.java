package edu.unc.mapseq.pipeline.rnaseq;

import java.util.Timer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class RNASeqPipelineExecutorService {

    private final Logger logger = LoggerFactory.getLogger(RNASeqPipelineExecutorService.class);

    private final Timer mainTimer = new Timer();

    private RNASeqPipelineExecutorTask task;

    private Long period = Long.valueOf(5);

    public RNASeqPipelineExecutorService() {
        super();
    }

    public void start() throws Exception {
        logger.info("ENTERING start()");
        long delay = 1 * 60 * 1000;
        mainTimer.scheduleAtFixedRate(task, delay, period * 60 * 1000);
    }

    public void stop() throws Exception {
        logger.info("ENTERING stop()");
        mainTimer.purge();
        mainTimer.cancel();
    }

    public RNASeqPipelineExecutorTask getTask() {
        return task;
    }

    public void setTask(RNASeqPipelineExecutorTask task) {
        this.task = task;
    }

    public Long getPeriod() {
        return period;
    }

    public void setPeriod(Long period) {
        this.period = period;
    }

}
