package edu.unc.mapseq.messaging.rnaseq;

import javax.jms.Connection;
import javax.jms.ConnectionFactory;
import javax.jms.Destination;
import javax.jms.MessageConsumer;
import javax.jms.Session;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class RNASeqMessageService {

    private final Logger logger = LoggerFactory.getLogger(RNASeqMessageService.class);

    private Connection connection;

    private Session session;

    private ConnectionFactory connectionFactory;

    private RNASeqMessageListener messageListener;

    private String destinationName;

    public void start() throws Exception {
        logger.info("ENTERING start()");
        this.connection = connectionFactory.createConnection();
        this.session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
        Destination destination = this.session.createQueue(this.destinationName);
        MessageConsumer consumer = session.createConsumer(destination);
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

    public RNASeqMessageListener getMessageListener() {
        return messageListener;
    }

    public void setMessageListener(RNASeqMessageListener messageListener) {
        this.messageListener = messageListener;
    }

    public String getDestinationName() {
        return destinationName;
    }

    public void setDestinationName(String destinationName) {
        this.destinationName = destinationName;
    }

}
