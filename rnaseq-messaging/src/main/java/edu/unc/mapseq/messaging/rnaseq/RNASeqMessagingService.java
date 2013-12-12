package edu.unc.mapseq.messaging.rnaseq;

import javax.jms.Connection;
import javax.jms.ConnectionFactory;
import javax.jms.Destination;
import javax.jms.MessageConsumer;
import javax.jms.Session;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.workflow.WorkflowBeanService;

public class RNASeqMessagingService {

    private final Logger logger = LoggerFactory.getLogger(RNASeqMessagingService.class);

    private Connection connection;

    private Session session;

    private ConnectionFactory connectionFactory;

    private WorkflowBeanService workflowBeanService;

    private String destinationName;

    public void start() throws Exception {
        logger.info("ENTERING start()");
        this.connection = connectionFactory.createConnection();
        this.session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
        Destination destination = this.session.createQueue(this.destinationName);
        MessageConsumer consumer = session.createConsumer(destination);
        RNASeqMessageListener messageListener = new RNASeqMessageListener();
        messageListener.setWorkflowBeanService(workflowBeanService);
        consumer.setMessageListener(messageListener);
        this.connection.start();
    }

    public void stop() throws Exception {
        logger.info("ENTERING stop()");
        if (this.session != null) {
            this.session.close();
        }
        if (this.connection != null) {
            this.connection.stop();
            this.connection.close();
        }
    }

    public ConnectionFactory getConnectionFactory() {
        return connectionFactory;
    }

    public void setConnectionFactory(ConnectionFactory connectionFactory) {
        this.connectionFactory = connectionFactory;
    }

    public WorkflowBeanService getWorkflowBeanService() {
        return workflowBeanService;
    }

    public void setWorkflowBeanService(WorkflowBeanService workflowBeanService) {
        this.workflowBeanService = workflowBeanService;
    }

    public String getDestinationName() {
        return destinationName;
    }

    public void setDestinationName(String destinationName) {
        this.destinationName = destinationName;
    }

}
