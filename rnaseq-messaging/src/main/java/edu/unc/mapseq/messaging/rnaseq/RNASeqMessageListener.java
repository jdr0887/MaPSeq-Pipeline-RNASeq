package edu.unc.mapseq.messaging.rnaseq;

import java.util.HashSet;
import java.util.Set;

import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageListener;
import javax.jms.TextMessage;

import org.apache.commons.lang.StringUtils;
import org.codehaus.jettison.json.JSONArray;
import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Account;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.dao.model.WorkflowPlan;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunStatusType;
import edu.unc.mapseq.pipeline.EntityUtil;
import edu.unc.mapseq.pipeline.rnaseq.RNASeqPipelineBeanService;

public class RNASeqMessageListener implements MessageListener {

    private final Logger logger = LoggerFactory.getLogger(RNASeqMessageListener.class);

    private RNASeqPipelineBeanService pipelineBeanService;

    public RNASeqMessageListener() {
        super();
    }

    @Override
    public void onMessage(Message message) {
        logger.debug("ENTERING onMessage(Message)");

        String messageValue = null;

        try {
            if (message instanceof TextMessage) {
                logger.debug("received TextMessage");
                TextMessage textMessage = (TextMessage) message;
                messageValue = textMessage.getText();
            }
        } catch (JMSException e2) {
            e2.printStackTrace();
        }

        if (StringUtils.isEmpty(messageValue)) {
            logger.warn("message value is empty");
            return;
        }

        logger.info("messageValue: {}", messageValue);

        JSONObject jsonMessage = null;

        try {
            jsonMessage = new JSONObject(messageValue);
            if (!jsonMessage.has("entities") || !jsonMessage.has("account_name")) {
                logger.error("json lacks entities or account_name");
                return;
            }
        } catch (JSONException e) {
            logger.error("BAD JSON format", e);
            return;
        }

        SequencerRun sequencerRun = null;
        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();
        WorkflowRun workflowRun = null;
        Account account = null;

        try {

            String accountName = jsonMessage.getString("account_name");

            try {
                account = this.pipelineBeanService.getMaPSeqDAOBean().getAccountDAO().findByName(accountName);
            } catch (MaPSeqDAOException e) {
            }

            if (account == null) {
                logger.error("Must register account first");
                return;
            }

            JSONArray entityArray = jsonMessage.getJSONArray("entities");

            for (int i = 0; i < entityArray.length(); ++i) {

                JSONObject entityJSONObject = entityArray.getJSONObject(i);

                if (entityJSONObject.has("entity_type")) {

                    String entityType = entityJSONObject.getString("entity_type");

                    if ("Sequencer run".equals(entityType) || SequencerRun.class.getSimpleName().equals(entityType)) {
                        sequencerRun = EntityUtil.getSequencerRun(pipelineBeanService.getMaPSeqDAOBean(),
                                entityJSONObject);
                    }

                    if ("HTSF Sample".equals(entityType) || HTSFSample.class.getSimpleName().equals(entityType)) {
                        HTSFSample htsfSample = EntityUtil.getHTSFSample(pipelineBeanService.getMaPSeqDAOBean(),
                                entityJSONObject);
                        htsfSampleSet.add(htsfSample);
                    }

                    if ("Workflow run".equals(entityType) || WorkflowRun.class.getSimpleName().equals(entityType)) {
                        workflowRun = EntityUtil.getWorkflowRun(pipelineBeanService.getMaPSeqDAOBean(), "RNASeq",
                                entityJSONObject, account);
                    }

                }

            }
        } catch (JSONException e1) {
            e1.printStackTrace();
            return;
        }

        if (workflowRun == null) {
            logger.warn("Invalid JSON...not running anything");
            return;
        }

        if (sequencerRun == null && htsfSampleSet.size() == 0) {
            logger.warn("Invalid JSON...not running anything");
            workflowRun.setStatus(WorkflowRunStatusType.FAILED);
        }

        try {
            Long workflowRunId = pipelineBeanService.getMaPSeqDAOBean().getWorkflowRunDAO().save(workflowRun);
            workflowRun.setId(workflowRunId);
        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

        try {
            WorkflowPlan workflowPlan = new WorkflowPlan();
            workflowPlan.setWorkflowRun(workflowRun);
            if (htsfSampleSet.size() > 0) {
                workflowPlan.setHTSFSamples(htsfSampleSet);
            }
            if (sequencerRun != null) {
                workflowPlan.setSequencerRun(sequencerRun);
            }
            this.pipelineBeanService.getMaPSeqDAOBean().getWorkflowPlanDAO().save(workflowPlan);
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
