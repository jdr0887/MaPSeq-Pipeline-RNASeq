package edu.unc.mapseq.executor.rnaseq;

import java.util.Timer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class RNASeqWorkflowExecutorService {

    private final Logger logger = LoggerFactory.getLogger(RNASeqWorkflowExecutorService.class);

    private final Timer mainTimer = new Timer();

    private RNASeqWorkflowExecutorTask task;

    private Long period = Long.valueOf(5);

    public RNASeqWorkflowExecutorService() {
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

    public RNASeqWorkflowExecutorTask getTask() {
        return task;
    }

    public void setTask(RNASeqWorkflowExecutorTask task) {
        this.task = task;
    }

    public Long getPeriod() {
        return period;
    }

    public void setPeriod(Long period) {
        this.period = period;
    }

}
