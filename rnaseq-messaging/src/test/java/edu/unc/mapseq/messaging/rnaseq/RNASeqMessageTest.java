package edu.unc.mapseq.messaging.rnaseq;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;
import org.junit.Test;

public class RNASeqMessageTest {

    @Test
    public void testQueue() {

        ActiveMQConnectionFactory connectionFactory = new ActiveMQConnectionFactory(String.format("nio://%s:61616",
                "gnet641.its.unc.edu"));

        Connection connection = null;
        Session session = null;
        try {
            connection = connectionFactory.createConnection();
            session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("queue/rnaseq");
            MessageProducer producer = session.createProducer(destination);
            producer.setDeliveryMode(DeliveryMode.PERSISTENT);
            String format = "{\"entities\":[{\"entityType\":\"Sample\",\"id\":\"%1$d\"},{\"entityType\":\"WorkflowRun\",\"name\":\"RNASeq-CALGB-%1$d\"}]}";
            // producer.send(session.createTextMessage(String.format(format, 397056L)));

            Integer[] sampleIdentifierArray = new Integer[] { 397057, 397058, 397059, 397060, 397061, 397062, 397063,
                    397064, 397065, 397066, 397067, 397068, 397069, 397070, 397071, 397072, 397073, 397074, 397075,
                    397076, 397077, 397078, 397079, 397172, 397173, 397174, 397176, 397177, 397179, 397180, 397182,
                    397183, 397184, 397186, 397187, 397189, 397190, 397192, 397193, 397194, 397196, 397197, 397199,
                    397201, 397203, 397204, 397205 };

            for (Integer sampleId : sampleIdentifierArray) {
                producer.send(session.createTextMessage(String.format(format, sampleId)));
            }

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
