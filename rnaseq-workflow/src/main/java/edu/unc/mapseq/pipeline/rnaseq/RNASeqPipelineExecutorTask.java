package edu.unc.mapseq.pipeline.rnaseq;

import java.util.Date;
import java.util.List;
import java.util.TimerTask;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.WorkflowDAO;
import edu.unc.mapseq.dao.WorkflowPlanDAO;
import edu.unc.mapseq.dao.WorkflowRunDAO;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowPlan;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.pipeline.PipelineExecutor;

public class RNASeqPipelineExecutorTask extends TimerTask {

    private final Logger logger = LoggerFactory.getLogger(RNASeqPipelineExecutorTask.class);

    private static final RNASeqPipelineTPE tpe = new RNASeqPipelineTPE();

    private RNASeqPipelineBeanService pipelineBeanService;

    public RNASeqPipelineExecutorTask() {
        super();
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");
        WorkflowDAO workflowDAO = this.pipelineBeanService.getMaPSeqDAOBean().getWorkflowDAO();
        WorkflowRunDAO workflowRunDAO = this.pipelineBeanService.getMaPSeqDAOBean().getWorkflowRunDAO();
        WorkflowPlanDAO workflowPlanDAO = this.pipelineBeanService.getMaPSeqDAOBean().getWorkflowPlanDAO();

        try {
            Workflow workflow = workflowDAO.findByName("RNASeq");
            List<WorkflowPlan> workflowPlanList = workflowPlanDAO.findEnqueued(workflow.getId());

            if (workflowPlanList != null && workflowPlanList.size() > 0) {

                logger.info("dequeuing {} WorkflowPlans", workflowPlanList.size());
                for (WorkflowPlan workflowPlan : workflowPlanList) {

                    RNASeqPipeline pipeline = new RNASeqPipeline();

                    WorkflowRun workflowRun = workflowPlan.getWorkflowRun();
                    workflowRun.setVersion(pipeline.getVersion());
                    workflowRun.setDequeuedDate(new Date());
                    workflowRunDAO.save(workflowRun);

                    // not sure if this is really necessary
                    workflowPlan.setWorkflowRun(workflowRun);

                    pipeline.setPipelineBeanService(pipelineBeanService);
                    pipeline.setWorkflowPlan(workflowPlan);
                    tpe.submit(new PipelineExecutor(pipeline));

                }

            }

        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

    }

    public RNASeqPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(RNASeqPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

}
