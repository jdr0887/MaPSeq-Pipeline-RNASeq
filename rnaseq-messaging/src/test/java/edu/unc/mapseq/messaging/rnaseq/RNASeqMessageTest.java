package edu.unc.mapseq.messaging.rnaseq;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MessageProducer;
import javax.jms.Session;
import javax.jms.TextMessage;

import org.apache.activemq.ActiveMQConnectionFactory;
import org.junit.Test;

public class RNASeqMessageTest {

    @Test
    public void testQueue() {
        // String host = "localhost";
        String host = "biodev1.its.unc.edu";
        ActiveMQConnectionFactory connectionFactory = new ActiveMQConnectionFactory("tcp://" + host + ":61616");

        Connection connection = null;
        Session session = null;
        try {
            connection = connectionFactory.createConnection();
            session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("queue/RNASeq");
            MessageProducer producer = session.createProducer(destination);
            producer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);
            // String value =
            // "{\"account_name\":\"jreilly\",\"entities\":[{\"entity_type\":\"Sequencer run\",\"guid\":\"45900\"},{\"entity_type\":\"Workflow run\",\"name\":\"RNASeq-jreilly-test1\"}]}";
            String value = "{\"account_name\":\"jreilly\",\"entities\":[{\"entity_type\":\"HTSF Sample\",\"guid\":\"45906\"},{\"entity_type\":\"Workflow run\",\"name\":\"RNASeq-jreilly-test1\"}]}";
            TextMessage message = session.createTextMessage(value);
            producer.send(message);
        } catch (JMSException e) {
            e.printStackTrace();
        } finally {
            try {
                session.close();
                connection.close();
            } catch (JMSException e) {
                e.printStackTrace();
            }
        }

    }

}
